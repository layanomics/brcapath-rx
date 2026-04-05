/* ============================================================
   BRCAPath-Rx  |  Bayes Demo  |  script.js
   Loads data.json, drives the interactive demo.
   ============================================================ */

var DATA = null;  // var so inline <script> tags can read it via window.DATA

/* ---------- Bootstrap ---------- */
fetch('data.json')
  .then(r => r.json())
  .then(d => {
    DATA = d;
    initDemo();
  })
  .catch(() => {
    document.getElementById('demo-error').style.display = 'block';
  });

/* ---------- State ---------- */
const state = { er: true, pr: true, her2: true, grade: true };

function lookupKey() {
  return `${+state.er}${+state.pr}${+state.her2}${+state.grade}`;
}

/* ---------- Init ---------- */
function initDemo() {
  // Wire toggle buttons
  document.querySelectorAll('.toggle-group button').forEach(btn => {
    btn.addEventListener('click', () => {
      const field = btn.dataset.field;
      const val   = btn.dataset.val === 'true';
      state[field] = val;
      // Update button active states
      btn.closest('.toggle-group').querySelectorAll('button').forEach(b => {
        b.classList.toggle('active', b.dataset.val === btn.dataset.val);
      });
      render();
    });
  });
  render();
}

/* ---------- Render ---------- */
function render() {
  if (!DATA) return;

  const key = lookupKey();
  const posteriors = DATA.posteriors[key];
  const priors     = DATA.priors;
  const likes      = DATA.likelihoods;

  // --- Compute combined likelihood ratio for display ---
  // P(evidence | LumA) — naive Bayes product
  const L = likes['LumA'];
  const lukVal =
    (state.er    ? L.er_pos   : 1 - L.er_pos)  *
    (state.pr    ? L.pr_pos   : 1 - L.pr_pos)  *
    (state.her2  ? L.her2_neg : 1 - L.her2_neg)*
    (state.grade ? L.low_grade: 1 - L.low_grade);

  // Prior, combined likelihood, posterior for LumA
  const prior     = priors['LumA'];
  const posterior = posteriors['LumA'];

  // Display pills
  setInner('val-prior',    pct(prior));
  setInner('val-prior-raw', `${DATA.counts['LumA']} / ${DATA.n_subtype_total}`);
  setInner('val-like',     pct(lukVal));
  setInner('val-post',     pct(posterior));

  // Directional arrow
  const arrow = posterior > prior + 0.02 ? '↑' : posterior < prior - 0.02 ? '↓' : '≈';
  const arrowEl = document.getElementById('post-arrow');
  if (arrowEl) {
    arrowEl.textContent = arrow;
    arrowEl.className = 'pill-sub';
    arrowEl.style.color = posterior > prior + 0.02 ? '#059669' : posterior < prior - 0.02 ? '#dc2626' : '#6b7280';
  }

  // Interpretation sentence
  const interp = buildInterpretation(prior, posterior, state);
  setInner('interp-text', interp);

  // Bar chart
  renderBars(posteriors);
}

function renderBars(posteriors) {
  const container = document.getElementById('bar-chart');
  if (!container) return;

  const subtypeClasses = {
    'LumA': 'luma', 'LumB': 'lumb', 'Her2': 'her2',
    'Basal': 'basal', 'Normal': 'normal', 'claudin-low': 'claudin'
  };
  const subtypeLabels = DATA.subtype_labels;

  container.innerHTML = '';
  const sorted = [...DATA.subtypes].sort((a, b) => posteriors[b] - posteriors[a]);

  sorted.forEach(s => {
    const p = posteriors[s];
    const row = document.createElement('div');
    row.className = 'bar-row' + (s === 'LumA' ? ' highlight' : '');

    const name = document.createElement('div');
    name.className = 'bar-name';
    name.textContent = subtypeLabels[s] || s;

    const track = document.createElement('div');
    track.className = 'bar-track';
    const fill = document.createElement('div');
    fill.className = `bar-fill ${subtypeClasses[s] || ''}`;
    fill.style.width = `${Math.round(p * 100)}%`;
    track.appendChild(fill);

    const pctEl = document.createElement('div');
    pctEl.className = 'bar-pct';
    pctEl.textContent = pct(p);

    row.appendChild(name);
    row.appendChild(track);
    row.appendChild(pctEl);
    container.appendChild(row);
  });
}

/* ---------- Helpers ---------- */
function pct(v) {
  return (v * 100).toFixed(1) + '%';
}
function setInner(id, html) {
  const el = document.getElementById(id);
  if (el) el.textContent = html;
}

function buildInterpretation(prior, posterior, s) {
  const erStr    = s.er    ? 'ER+' : 'ER−';
  const prStr    = s.pr    ? 'PR+' : 'PR−';
  const her2Str  = s.her2  ? 'HER2−' : 'HER2+';
  const gradeStr = s.grade ? 'low grade' : 'high grade';
  const change   = posterior > prior + 0.02 ? 'increases to'
                 : posterior < prior - 0.02 ? 'drops to'
                 : 'remains near';
  return `Given ${erStr}, ${prStr}, ${her2Str}, ${gradeStr}: the probability of Luminal A ${change} ${pct(posterior)} (from a prior of ${pct(prior)}).`;
}
