// Recreates the Array-CGH quick guide's screenshots (src/assets/cgh-guide/step*.webp) from the demo inputs here.
// Synthetic data only. Run from the repo root after `npm run build`: node docs/cgh-demo/make-guide-shots.cjs
// Needs playwright (chromium) and the network (UCSC, PanelApp, gnomAD annotate the calls).
const { chromium } = require('playwright');
const fs = require('fs'), path = require('path');
const D = path.resolve('docs/cgh-demo'), OUT = path.resolve('src/assets/cgh-guide');
const PAGE = 'file://' + path.resolve('sashimi-viewer.html') + '#sv=chr1';
const MAX_W = 900, Q = 0.7;

(async () => {
  fs.mkdirSync(OUT, { recursive: true });
  const b = await chromium.launch();
  const p = await b.newPage({ viewport: { width: 1400, height: 900 }, deviceScaleFactor: 1 });
  p.on('pageerror', e => console.error('page error:', e.message));
  await p.goto(PAGE);
  await p.click('button:has-text("Array-CGH")');
  // the red outline on the controls a step is about
  await p.addStyleTag({ content: '[data-guide-mark]{outline:3px solid #dc2626!important;outline-offset:2px;border-radius:8px}' });
  const mark = async sel => p.locator(sel).evaluateAll(es => es.forEach(e => e.setAttribute('data-guide-mark', '')));
  const unmark = () => p.evaluate(() => document.querySelectorAll('[data-guide-mark]').forEach(e => e.removeAttribute('data-guide-mark')));
  const box = sel => p.locator(sel).first().boundingBox();
  // a clipped PNG, written as WebP (scaled to MAX_W at most) by the page's own canvas
  const shot = async (name, clip) => {
    const png = (await p.screenshot({ clip })).toString('base64');
    const url = await p.evaluate(async ([src, maxW, q]) => {
      const img = new Image(); img.src = 'data:image/png;base64,' + src; await img.decode();
      const k = Math.min(1, maxW / img.width), c = document.createElement('canvas');
      c.width = Math.round(img.width * k); c.height = Math.round(img.height * k);
      const g = c.getContext('2d'); g.imageSmoothingQuality = 'high'; g.drawImage(img, 0, 0, c.width, c.height);
      return c.toDataURL('image/webp', q);
    }, [png, MAX_W, Q]);
    const buf = Buffer.from(url.split(',')[1], 'base64');
    fs.writeFileSync(path.join(OUT, name + '.webp'), buf);
    console.log(name, buf.length, 'bytes');
    await unmark();
  };
  const input = '[data-cgh-page] input[data-cgh-input] >> nth=0';

  // 1. the empty page: Add samples, or drop files
  await mark('[data-cgh-page] label:has(input[data-cgh-input])');
  await shot('step1', { x: 0, y: 0, width: 1400, height: 310 });

  // 2. the import question
  await p.setInputFiles(input, [path.join(D, 'DEMO-01.txt'), path.join(D, 'DEMO-02.txt')]);
  const dlg = '[role=dialog][aria-label="Minimum number of probes"] > div';
  await p.waitForSelector(dlg);
  await mark(`${dlg} span:has(input[aria-label="Minimum size (slider)"])`);
  await mark(`${dlg} span:has(> input[aria-label="Minimum number of probes"])`);
  await mark(`${dlg} table`);
  const d = await box(dlg);
  await shot('step2', { x: d.x - 8, y: d.y - 8, width: d.width + 16, height: d.height + 16 });
  await p.click('button:has-text("Analyse")');
  // the annotations: no cell still waiting (…)
  await p.waitForFunction(() => { const r = [...document.querySelectorAll('[data-cgh-row] td')]; return r.length && !r.some(td => /^…$|…$/.test(td.textContent.trim())); }, null, { timeout: 180000, polling: 1000 });
  await p.waitForTimeout(2000);

  // 5. the lab benign list
  await p.setInputFiles('input[data-lab-input]', path.join(D, 'lab-benign-demo.bed'));
  await p.click('button[role=tab]:has-text("All patients")');
  await p.waitForTimeout(800);
  await mark('label:has(input[data-lab-input])'); await mark('[data-lab-chip]');
  await mark('[data-cgh-row] span:text-is("lab benign")'); await mark('button:has-text("→ Benign?")');
  const fr = await box('[data-cgh-frame]');
  await shot('step5', { x: fr.x, y: fr.y, width: fr.width, height: 400 });

  // 4. the filters, More filters open
  await p.click('button:has-text("More filters")');
  await p.waitForTimeout(400);
  await mark('label:has-text("series")'); await mark('label:has-text("lab list")');
  await mark('input[aria-label="Filter the calls by gene, band or chromosome"]');
  await shot('step4', { x: fr.x, y: fr.y + 38, width: fr.width, height: 196 });
  await p.mouse.click(5, 880);

  // 3. the table: a patient's tab, a call clicked, its detail under it
  await p.setViewportSize({ width: 1400, height: 1300 });
  await p.click('button:has-text("DEMO-02")');
  await p.click('[data-cgh-row]:has-text("1q21.1")');
  await p.waitForTimeout(5000);
  await mark('[role=tablist]');
  const f3 = await box('[data-cgh-frame]');
  await shot('step3', { x: f3.x, y: f3.y, width: f3.width, height: Math.min(f3.height, 640) });

  // 3b. the browser below the table: the call clicked, drawn with the SV browser's tracks
  // the whole card (its ISCN line, then the browser); ClinVar SV off for the shot, so the reference and the rebuilt
  // chromosome fit in the image
  const card = p.locator('[data-cgh-page] .rounded-xl:has(span.font-mono:has-text("arr["))').last();
  const cvSv = card.locator('button:has-text("ClinVar SV")').first();
  await cvSv.click(); await p.waitForTimeout(800);
  // the card's top just under the page's sticky bar (48 px), so its ISCN line and the browser's header show
  await card.evaluate(el => { const pg = document.querySelector('[data-cgh-page]'); pg.scrollTop += el.getBoundingClientRect().top - 56; });
  await p.waitForTimeout(4000);
  const fb = await card.boundingBox();
  await mark('[data-cgh-page] span.font-mono:has-text("arr[")');
  await shot('step3b', { x: fb.x, y: fb.y - 4, width: fb.width, height: Math.min(fb.height + 4, 780) });
  await cvSv.click(); await p.waitForTimeout(300);
  await p.evaluate(() => document.querySelector('[data-cgh-page]').scrollTo(0, 0));

  // 6. classify, note, export
  const row = '[data-cgh-row]:has-text("1q21.1")';
  await p.selectOption(`${row} select`, 'Uncertain significance');
  await p.click(`${row} button[aria-label="Note"]`);
  await p.waitForTimeout(300);
  await p.keyboard.type('Demo note: parents to test');
  console.log('note:', await p.inputValue('textarea[aria-label="Note on this call"]'));
  await p.click('th:has-text("Dosage")'); // the note kept (blur)
  await p.waitForTimeout(500);
  await mark(`${row} select`); await mark(`${row} button[aria-label="Note"]`); await mark('button:has-text("Export classified table")');
  const f6 = await box('[data-cgh-frame]');
  await shot('step6', { x: f6.x + f6.width - 560, y: f6.y, width: 560, height: 330 });
  await b.close();
})();
