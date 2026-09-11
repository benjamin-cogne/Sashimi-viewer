// Archives the repository's GitHub traffic beyond the 14-day window GitHub keeps.
// Run by .github/workflows/traffic.yml weekly; usable by hand: GH_TOKEN=… node scripts/archive-traffic.mjs
//
// Writes, under docs/stats/ (created if needed):
//   views.csv      date, views, unique_visitors          one row per day (rows of the last 14 days are refreshed)
//   clones.csv     date, clones, unique_cloners          idem
//   referrers.csv  collected, referrer, views, unique_visitors     snapshot of the top referrers at each run
//   paths.csv      collected, path, title, views, unique_visitors  snapshot of the most visited repository pages
//   snapshot.csv   collected, stars, forks, watchers, release_downloads   one row per run
// Traffic endpoints need push access (the workflow token has it). Public endpoints still work without a token.
import { readFileSync, writeFileSync, mkdirSync, existsSync } from 'node:fs';
import { join } from 'node:path';

const REPO = process.env.GITHUB_REPOSITORY || 'benjamin-cogne/Sashimi-viewer';
const TOKEN = process.env.GH_TOKEN || process.env.GITHUB_TOKEN || '';
const OUT = process.env.STATS_DIR || 'docs/stats';
const API = 'https://api.github.com';
const today = new Date().toISOString().slice(0, 10);

async function get(path) {
  const headers = { Accept: 'application/vnd.github+json', 'User-Agent': 'sashimi-viewer-traffic-archive', 'X-GitHub-Api-Version': '2022-11-28' };
  if (TOKEN) headers.Authorization = `Bearer ${TOKEN}`;
  const r = await fetch(`${API}${path}`, { headers });
  if (!r.ok) { console.warn(`skip ${path}: HTTP ${r.status}`); return null; }
  return r.json();
}

const csvEscape = v => { const s = String(v ?? ''); return /[",\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s; };
/** Rows of a CSV written by writeCsv (header dropped); handles quoted cells with commas and doubled quotes. */
function readCsv(file) {
  if (!existsSync(file)) return [];
  const rows = [];
  for (const line of readFileSync(file, 'utf8').split('\n').slice(1)) {
    if (!line) continue;
    const cells = []; let cur = '', q = false;
    for (let i = 0; i < line.length; i++) {
      const c = line[i];
      if (q) { if (c === '"' && line[i + 1] === '"') { cur += '"'; i++; } else if (c === '"') q = false; else cur += c; }
      else if (c === '"') q = true;
      else if (c === ',') { cells.push(cur); cur = ''; }
      else cur += c;
    }
    cells.push(cur); rows.push(cells);
  }
  return rows;
}
function writeCsv(file, header, rows) {
  writeFileSync(file, [header.join(','), ...rows.map(r => r.map(csvEscape).join(','))].join('\n') + '\n');
}
/** Merge per-day rows: rows for dates present in the fresh data replace the archived ones (GitHub revises recent days). */
function mergeDaily(file, header, fresh) {
  const old = readCsv(file).filter(r => r.length >= header.length);
  const byDate = new Map(old.map(r => [r[0], r]));
  for (const r of fresh) byDate.set(r[0], r);
  const rows = [...byDate.values()].sort((a, b) => a[0].localeCompare(b[0]));
  writeCsv(file, header, rows);
  return rows.length;
}
function appendSnapshot(file, header, rows) {
  const old = readCsv(file).filter(r => r.length >= header.length && r[0] !== today); // rerunning the same day replaces that day
  writeCsv(file, header, [...old, ...rows]);
}

mkdirSync(OUT, { recursive: true });
const day = ts => String(ts).slice(0, 10);

const views = await get(`/repos/${REPO}/traffic/views`);
if (views) console.log(`views: ${mergeDaily(join(OUT, 'views.csv'), ['date', 'views', 'unique_visitors'], views.views.map(v => [day(v.timestamp), v.count, v.uniques]))} days archived`);
const clones = await get(`/repos/${REPO}/traffic/clones`);
if (clones) console.log(`clones: ${mergeDaily(join(OUT, 'clones.csv'), ['date', 'clones', 'unique_cloners'], clones.clones.map(v => [day(v.timestamp), v.count, v.uniques]))} days archived`);
const referrers = await get(`/repos/${REPO}/traffic/popular/referrers`);
if (referrers) { appendSnapshot(join(OUT, 'referrers.csv'), ['collected', 'referrer', 'views', 'unique_visitors'], referrers.map(r => [today, r.referrer, r.count, r.uniques])); console.log(`referrers: ${referrers.length}`); }
const paths = await get(`/repos/${REPO}/traffic/popular/paths`);
if (paths) { appendSnapshot(join(OUT, 'paths.csv'), ['collected', 'path', 'title', 'views', 'unique_visitors'], paths.map(p => [today, p.path, p.title, p.count, p.uniques])); console.log(`paths: ${paths.length}`); }

const repo = await get(`/repos/${REPO}`);
const releases = (await get(`/repos/${REPO}/releases?per_page=100`)) || [];
const downloads = releases.reduce((n, rel) => n + (rel.assets || []).reduce((m, a) => m + (a.download_count || 0), 0), 0);
if (repo) {
  appendSnapshot(join(OUT, 'snapshot.csv'), ['collected', 'stars', 'forks', 'watchers', 'release_downloads'], [[today, repo.stargazers_count, repo.forks_count, repo.subscribers_count, downloads]]);
  console.log(`snapshot: ${repo.stargazers_count} stars, ${repo.forks_count} forks, ${repo.subscribers_count} watchers, ${downloads} release downloads`);
}
