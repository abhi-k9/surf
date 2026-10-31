#!/usr/bin/env python3
"""Loopback-only UI bridge. Geometry is computed exclusively by persistent C++ workers."""
import argparse
import atexit
import json
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import signal
import subprocess
import threading
import time
from urllib.parse import urlparse, parse_qs

ROOT = Path(__file__).resolve().parents[1]
WORKERS = {}
LOCKS = {5: threading.Lock(), 9: threading.Lock()}
CLI = None

def worker(n, request):
    with LOCKS[n]:
        process = WORKERS.get(n)
        if process is None or process.poll() is not None:
            process = subprocess.Popen([str(CLI), str(ROOT / 'handoff'), str(n), 'worker'],
                                       stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=subprocess.DEVNULL, text=True, encoding='utf-8')
            WORKERS[n] = process
        process.stdin.write(json.dumps(request, allow_nan=False, separators=(',', ':')) + '\n')
        process.stdin.flush()
        reply = process.stdout.readline(32 * 1024 * 1024)
        if not reply.endswith('\n'):
            process.kill()
            raise ValueError('C++ worker failed or response exceeded limit')
        result = json.loads(reply)
        if 'error' in result:
            raise ValueError(result['error'])
        return result

@atexit.register
def close_workers():
    for p in WORKERS.values():
        if p.poll() is None:
            p.terminate()
            try:
                p.wait(timeout=3)
            except subprocess.TimeoutExpired:
                p.kill()
                p.wait()

class Handler(BaseHTTPRequestHandler):
    def setup(self):
        super().setup()
        self.connection.settimeout(15)

    def log_message(self, *_):
        pass

    def send(self, code, body, content='application/json'):
        if isinstance(body, (dict, list)):
            body = json.dumps(body, allow_nan=False, separators=(',', ':')).encode()
        self.send_response(code)
        self.send_header('Content-Type', content)
        self.send_header('Content-Length', str(len(body)))
        self.send_header('Cache-Control', 'no-store')
        self.send_header('X-Content-Type-Options', 'nosniff')
        self.send_header('Content-Security-Policy', "default-src 'self'; script-src 'self'; style-src 'self'; img-src 'self' data:; connect-src 'self'; frame-ancestors 'none'; base-uri 'none'")
        self.end_headers()
        self.wfile.write(body)

    def allowed(self):
        # Mitigate DNS rebinding and cross-origin requests to the local bridge.
        port = self.server.server_port
        allowed = {f'127.0.0.1:{port}', f'localhost:{port}'}
        host = self.headers.get('Host', '')
        origin = self.headers.get('Origin')
        return host in allowed and (origin is None or origin == f'http://{host}')

    def do_GET(self):
        if not self.allowed():
            return self.send(403, {'error': 'local same-origin requests required'})
        url = urlparse(self.path)
        if url.path in ('/', '/app.js', '/style.css'):
            filename = {'/': 'index.html', '/app.js': 'app.js', '/style.css': 'style.css'}[url.path]
            return self.send(200, (ROOT / 'viewer' / filename).read_bytes(),
                             {'/': 'text/html; charset=utf-8', '/app.js': 'text/javascript; charset=utf-8', '/style.css': 'text/css; charset=utf-8'}[url.path])
        if url.path == '/api/meta':
            try:
                n = int(parse_qs(url.query).get('valence', ['5'])[0])
                if n not in (5, 9):
                    raise ValueError('unsupported fan')
                schema = json.loads((ROOT / 'handoff/fixtures' / f'fan{n}_schema.json').read_text())
                scenes = {s: json.loads((ROOT / 'handoff/scenes' / f'fan{n}_{s}.json').read_text()) for s in ('plane', 'rational')}
                inventory = json.loads(subprocess.check_output([str(CLI), str(ROOT / 'handoff'), str(n), 'inventory'], timeout=10))
                return self.send(200, {'schema': schema, 'scenes': scenes, 'supports': inventory['supports']})
            except (ValueError, KeyError) as e:
                return self.send(400, {'error': str(e)})
        self.send(404, {'error': 'not found'})

    def do_POST(self):
        if not self.allowed():
            return self.send(403, {'error': 'local same-origin requests required'})
        if self.path != '/api/mesh':
            return self.send(404, {'error': 'not found'})
        try:
            if self.headers.get('Content-Type') != 'application/json':
                raise ValueError('application/json required')
            length = int(self.headers.get('Content-Length', '-1'))
            if not 0 < length <= 2 * 1024 * 1024:
                raise ValueError('body must be 1 byte to 2 MiB')
            request = json.loads(self.rfile.read(length), parse_constant=lambda _: (_ for _ in ()).throw(ValueError('nonfinite JSON')))
            n = request.get('valence')
            if type(n) is not int or n not in (5, 9):
                raise ValueError('unsupported fan')
            started = time.perf_counter()
            result = worker(n, request)
            result['bridge_ms'] = (time.perf_counter() - started) * 1000
            return self.send(200, result)
        except (ValueError, TypeError, KeyError, BrokenPipeError, RecursionError) as e:
            return self.send(400, {'error': str(e)})

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--cli', type=Path, default=ROOT / 'build/release/npms_cli')
    parser.add_argument('--port', type=int, default=8765)
    args = parser.parse_args()
    CLI = args.cli.resolve()
    if not CLI.is_file():
        parser.error('build npms_cli first, or supply --cli PATH')
    def stop_server(*_):
        raise KeyboardInterrupt
    signal.signal(signal.SIGTERM, stop_server)
    server = ThreadingHTTPServer(('127.0.0.1', args.port), Handler)
    server.daemon_threads = True
    print(f'NPMS local viewer on port {server.server_port}. Open 127.0.0.1:{server.server_port} in a browser. Ctrl+C stops.', flush=True)
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        pass
    finally:
        server.server_close()
        close_workers()
