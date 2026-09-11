# Usage statistics

Written by `.github/workflows/traffic.yml` (weekly, or on demand from the Actions tab) with
`scripts/archive-traffic.mjs`, because GitHub keeps only 14 days of traffic. Figures come from the
GitHub API and concern the repository on github.com, not the GitHub Pages site, which has no logs.

| File | Content |
|---|---|
| `views.csv` | repository page views and unique visitors per day |
| `clones.csv` | git clones and unique cloners per day |
| `referrers.csv` | top referring sites at each run (`collected` is the run date) |
| `paths.csv` | most visited repository pages at each run |
| `snapshot.csv` | stars, forks, watchers and cumulated release-asset downloads at each run |

Rows of the last 14 days are refreshed at each run, since GitHub revises recent counts.
