"""Browser interaction evidence; requires Playwright and Chromium, optional to core CI."""
import json
from pathlib import Path
import time
from playwright.sync_api import sync_playwright
ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'evidence'
OUT.mkdir(exist_ok=True)
with sync_playwright() as p:
    browser = p.chromium.launch(executable_path='/usr/bin/chromium', headless=True,
                                args=['--no-sandbox', '--use-angle=swiftshader', '--enable-unsafe-swiftshader'])
    page = browser.new_page(viewport={'width':1440,'height':1100}, device_scale_factor=1)
    errors = []
    page.on('pageerror', lambda e: errors.append(str(e)))
    page.goto('http://127.0.0.1:8765', wait_until='networkidle')
    def ready():
        page.wait_for_function('() => window.npmsState && window.npmsState().triangles > 0 && !window.npmsState().busy')
        assert not page.locator('#error').inner_text(), page.locator('#error').inner_text()
    def shot(name):
        page.screenshot(path=str(OUT / f'{name}.png'))
    ready()
    records=[]
    for n in ('5','9'):
        page.select_option('#fan', n)
        page.wait_for_function(f'() => document.getElementById("statFan").textContent==="0{n}" && !window.npmsState().busy')
        ready()
        page.uncheck('#cage')
        page.select_option('#preset','plane');ready()
        page.check('#ps');page.select_option('#mode','shade');shot(f'fan{n}-plane')
        page.select_option('#preset','rational');ready()
        before=page.evaluate('window.npmsState()')
        page.fill('#z',str(before['scene']['control_xyz'][0][2]+.5));page.click('#apply');ready()
        after=page.evaluate('window.npmsState()')
        assert after['scene']['control_xyz'][0][2]!=before['scene']['control_xyz'][0][2]
        page.click('#undo');ready();assert page.evaluate('window.npmsState().scene')==before['scene']
        page.click('#redo');ready();assert page.evaluate('window.npmsState().scene')==after['scene']
        page.fill('#weight','1.25');page.click('#apply');ready();assert page.evaluate('window.npmsState().scene.weights[0]')==1.25
        shot(f'fan{n}-edited')
        # Capture a parent group with deliberate duplicate sites and true influence.
        page.select_option('#parent','1');ready();page.select_option('#mode','influence');ready();shot(f'fan{n}-influence')
        page.select_option('#mode','mean');shot(f'fan{n}-curvature')
        # Group editing changes each matching child by the same amount.
        page.check('#group');before_group=page.evaluate('window.npmsState()')
        page.fill('#z',str(before_group['scene']['control_xyz'][before_group['selected']][2]+.1));page.click('#apply');ready()
        group_after=page.evaluate('window.npmsState()')
        assert sum(a!=b for a,b in zip(before_group['scene']['control_xyz'],group_after['scene']['control_xyz']))>1
        page.uncheck('#group')
        # Strict positive weight validation preserves last valid scene.
        previous=page.evaluate('window.npmsState().scene');page.fill('#weight','0');page.click('#apply');assert page.locator('#error').inner_text();assert page.evaluate('window.npmsState().scene')==previous
        page.fill('#weight','1');page.click('#apply');ready()
        # Save and reload exact scene, OBJ normals, camera navigation and picking.
        with page.expect_download() as download:
            page.click('#save')
        scene_file=OUT/f'fan{n}-saved-scene.json';download.value.save_as(scene_file)
        saved=json.loads(scene_file.read_text());assert saved==page.evaluate('window.npmsState().scene')
        page.click('#reset');ready();page.set_input_files('#file',str(scene_file));page.wait_for_function('(expected) => JSON.stringify(window.npmsState().scene.control_xyz) === JSON.stringify(expected.control_xyz) && JSON.stringify(window.npmsState().scene.weights) === JSON.stringify(expected.weights)',arg=saved);ready();assert saved['control_xyz']==page.evaluate('window.npmsState().scene.control_xyz');assert saved['weights']==page.evaluate('window.npmsState().scene.weights')
        with page.expect_download() as download:
            page.click('#obj')
        obj_file=OUT/f'fan{n}-edited.obj';download.value.save_as(obj_file);assert '\nvn ' in obj_file.read_text()
        page.mouse.move(900,400);page.mouse.down();page.mouse.move(1020,450,steps=20);page.mouse.up();page.mouse.wheel(0,100)
        page.click('#fitControls')
        pick=page.evaluate('''() => {const r=canvas.getBoundingClientRect(),a=axes(),p=scene.control_xyz[0].map((x,k)=>x-target[k]),z=distance-dot(p,a.forward);return {x:r.left+(1+1.9*(dot(p,a.right)+pan[0])/z/(r.width/r.height))*r.width/2,y:r.top+(1-1.9*(dot(p,a.up)+pan[1])/z)*r.height/2};}''')
        page.mouse.click(pick['x'],pick['y']);ready();assert page.evaluate('window.npmsState().selected')==0
        page.click('#fit')
        profile=page.evaluate('() => window.npmsProfileFrames(120)')
        records.append({'fan':n,'state':page.evaluate('window.npmsState()'),'profile':profile})
        page.click('#collapse');ready();assert page.evaluate('window.npmsState().invalid')>0;shot(f'fan{n}-invalid')
        page.click('#reset');ready()
    renderer=page.evaluate('''() => {const gl=document.getElementById('canvas').getContext('webgl2');const e=gl.getExtension('WEBGL_debug_renderer_info');return {version:gl.getParameter(gl.VERSION),renderer:e?gl.getParameter(e.UNMASKED_RENDERER_WEBGL):gl.getParameter(gl.RENDERER)};}''')
    assert not errors, errors
    (OUT/'viewer-tests.json').write_text(json.dumps({'passed':True,'browser':browser.version,'renderer':renderer,'records':records,'page_errors':errors},indent=2))
    print('PASS both fans: edit XYZ and weights, group edit, undo/redo, save/reload, OBJ, orbit/zoom/picking smoke, influence, valid/invalid curvature; 10 screenshots')
    browser.close()
