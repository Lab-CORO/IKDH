'use strict';

// ─── Module state ─────────────────────────────────────────────────────────────
let _frames      = null;   // last FK frames (7 × Float64Array 16)
let _Jspatial    = null;   // spatial (base-frame) Jacobian 6×6
let _transData   = null;   // { sigmas, vecs } from spatial J — used for 3D ellipsoid
let _rotData     = null;
let _showTransEl  = false;
let _showTransVec = false;
let _showRotEl    = false;
let _showRotVec   = false;
let _ellGroup    = null;   // THREE.Group added to scene

// ─── Init (call once after scene is ready) ───────────────────────────────────
function initAnalysis(threeScene) {
  _ellGroup = new THREE.Group();
  threeScene.add(_ellGroup);

  document.querySelector('.jac-eq-label').innerHTML =
    katex.renderToString('J =', {throwOnError: false});
}

function clearAnalysisState() {
  _frames = null; _Jspatial = null;
  _transData = null; _rotData = null;
  _showTransEl = false; _showTransVec = false;
  _showRotEl   = false; _showRotVec   = false;
  _rebuildScene();
}

// ─── Linear algebra ──────────────────────────────────────────────────────────
function cross3(a, b) {
  return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]];
}

// Evaluate the symbolic Jacobian rows from YAML at joint angles q (radians).
// Shorthand: s1=sin(q1), c1=cos(q1), ..., sin/cos are Math.sin/cos.
function evalJacobian(rows, q) {
  /* jshint evil: true */
  const [q1, q2, q3, q4, q5, q6] = q;
  const s1=Math.sin(q1), c1=Math.cos(q1);
  const s2=Math.sin(q2), c2=Math.cos(q2);
  const s3=Math.sin(q3), c3=Math.cos(q3);
  const s4=Math.sin(q4), c4=Math.cos(q4);
  const s5=Math.sin(q5), c5=Math.cos(q5);
  const s6=Math.sin(q6), c6=Math.cos(q6);
  const sin=Math.sin, cos=Math.cos;
  return rows.map(row => row.map(expr => eval(expr))); // eslint-disable-line no-eval
}

// A·Aᵀ for a 3×6 sub-matrix (3 rows, each 6 elements) → 3×3 flat row-major
function aat3(rows) {
  const A = new Array(9).fill(0);
  for (let r=0;r<3;r++) for (let c=0;c<3;c++) for (let k=0;k<6;k++) A[r*3+c]+=rows[r][k]*rows[c][k];
  return A;
}

// Jacobi eigendecomposition for 3×3 symmetric matrix. Returns {vals, vecs} descending.
function eig3sym(M) {
  let a=[...M], V=[1,0,0,0,1,0,0,0,1];
  for (let it=0;it<60;it++) {
    let p=0,q=1;
    if (Math.abs(a[2])>Math.abs(a[p*3+q])){p=0;q=2;}
    if (Math.abs(a[5])>Math.abs(a[p*3+q])){p=1;q=2;}
    if (Math.abs(a[p*3+q])<1e-14) break;
    const th=0.5*Math.atan2(2*a[p*3+q],a[q*3+q]-a[p*3+p]);
    const c=Math.cos(th),s=Math.sin(th),r=3-p-q;
    const na=[...a];
    na[p*3+p]=c*c*a[p*3+p]-2*s*c*a[p*3+q]+s*s*a[q*3+q];
    na[q*3+q]=s*s*a[p*3+p]+2*s*c*a[p*3+q]+c*c*a[q*3+q];
    na[p*3+q]=0;na[q*3+p]=0;
    na[r*3+p]=c*a[r*3+p]-s*a[r*3+q];na[p*3+r]=na[r*3+p];
    na[r*3+q]=s*a[r*3+p]+c*a[r*3+q];na[q*3+r]=na[r*3+q];
    na[r*3+r]=a[r*3+r]; a=na;
    const nV=[...V];
    for (let row=0;row<3;row++){
      nV[row*3+p]=c*V[row*3+p]-s*V[row*3+q];
      nV[row*3+q]=s*V[row*3+p]+c*V[row*3+q];
    }
    V=nV;
  }
  const ord=[0,1,2].sort((i,j)=>a[j*3+j]-a[i*3+i]);
  return { vals:ord.map(i=>Math.max(0,a[i*3+i])), vecs:ord.map(i=>[V[i],V[3+i],V[6+i]]) };
}

function _det6(J) {
  const A = J.map(r => [...r]);
  let s = 1;
  for (let c = 0; c < 6; c++) {
    let p = c;
    for (let r = c + 1; r < 6; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
    if (p !== c) { [A[c], A[p]] = [A[p], A[c]]; s = -s; }
    if (Math.abs(A[c][c]) < 1e-14) return 0;
    for (let r = c + 1; r < 6; r++) {
      const f = A[r][c] / A[c][c];
      for (let k = c; k < 6; k++) A[r][k] -= f * A[c][k];
    }
  }
  return s * A.reduce((d, r, i) => d * r[i], 1);
}

// ─── DOM helpers ─────────────────────────────────────────────────────────────
function _renderJacobianTable(J) {
  const rl = ['v<sub>x</sub>','v<sub>y</sub>','v<sub>z</sub>','ω<sub>x</sub>','ω<sub>y</sub>','ω<sub>z</sub>'];
  const cols = ['J1','J2','J3','J4','J5','J6'];
  document.getElementById('jac-col-hdr').innerHTML =
    cols.map(c => `<span>${c}</span>`).join('');
  document.getElementById('jac-row-labels').innerHTML =
    rl.map(r => `<span>${r}</span>`).join('');
  document.getElementById('jac-table').innerHTML =
    `<tbody>${J.map(row =>
      `<tr>${row.map(v=>`<td>${v.toFixed(3)}</td>`).join('')}</tr>`
    ).join('')}</tbody>`;
}

function _renderManipBlock(id, prefix, data) {
  const {sigmas} = data;
  const vol = sigmas[0]*sigmas[1]*sigmas[2];
  const isT = prefix === 't';
  const latexStr = isT
    ? '\\mu_t = \\sqrt{\\det(J_v J_v^T)} ='
    : '\\mu_o = \\sqrt{\\det(J_\\omega J_\\omega^T)} =';
  const formula = katex.renderToString(latexStr, {throwOnError: false});
  document.getElementById(id).innerHTML =
    `<div class="manip-row" style="font-size:0.82rem;align-items:baseline;gap:6px">
      ${formula}<span class="m-val">${vol.toExponential(3)}</span>
    </div>
    <div class="manip-row">
      <span class="m-lbl">σ₁</span><span class="m-val">${sigmas[0].toFixed(4)}</span>
      <span class="m-lbl">σ₂</span><span class="m-val">${sigmas[1].toFixed(4)}</span>
      <span class="m-lbl">σ₃</span><span class="m-val">${sigmas[2].toFixed(4)}</span>
      <label class="ell-check" style="margin-left:8px"><input type="checkbox" id="cb-vec-${prefix}"> Vectors</label>
      <label class="ell-check"><input type="checkbox" id="cb-ell-${prefix}"> Ellipsoid</label>
    </div>`;
  document.getElementById(`cb-ell-${prefix}`).addEventListener('change', e => {
    if (isT) _showTransEl  = e.target.checked;
    else     _showRotEl    = e.target.checked;
    _rebuildScene();
  });
  document.getElementById(`cb-vec-${prefix}`).addEventListener('change', e => {
    if (isT) _showTransVec = e.target.checked;
    else     _showRotVec   = e.target.checked;
    _rebuildScene();
  });
}

function _syncCheckboxes() {
  const te = document.getElementById('cb-ell-t');
  const tv = document.getElementById('cb-vec-t');
  const re = document.getElementById('cb-ell-r');
  const rv = document.getElementById('cb-vec-r');
  if (te) te.checked = _showTransEl;
  if (tv) tv.checked = _showTransVec;
  if (re) re.checked = _showRotEl;
  if (rv) rv.checked = _showRotVec;
}

// ─── Three.js objects ────────────────────────────────────────────────────────
function _makeEllipsoid(data, pos, color, visualRadius) {
  const {sigmas, vecs} = data;
  const maxS = Math.max(...sigmas, 1e-9);
  const s = sigmas.map(v => v/maxS * visualRadius);

  // Ensure right-handed basis (odd sort permutations give det = -1)
  const v0 = vecs[0], v1 = vecs[1], v2 = cross3(v0, v1);

  // Bake rotation+scale into geometry vertices to bypass quaternion extraction.
  const geo = new THREE.SphereGeometry(1, 24, 16);
  geo.applyMatrix4(new THREE.Matrix4().set(
    s[0]*v0[0], s[1]*v1[0], s[2]*v2[0], 0,
    s[0]*v0[1], s[1]*v1[1], s[2]*v2[1], 0,
    s[0]*v0[2], s[1]*v1[2], s[2]*v2[2], 0,
    0,          0,          0,          1
  ));

  const mesh = new THREE.Mesh(geo, new THREE.MeshPhongMaterial({
    color, opacity: 0.22, transparent: true,
    side: THREE.DoubleSide, depthWrite: false, shininess: 20
  }));
  mesh.position.copy(pos);
  return mesh;
}

function _makeAxisArrows(data, pos, color, ellRadius) {
  const {sigmas, vecs} = data;
  const maxS = Math.max(...sigmas, 1e-9);
  const dirs = [vecs[0], vecs[1], cross3(vecs[0], vecs[1])];
  const axLen = robotReach * 0.12;
  const headL = axLen * 0.28;
  const headW = axLen * 0.18;

  const group = new THREE.Group();
  dirs.forEach((dir, i) => {
    const len = sigmas[i] / maxS * ellRadius;
    if (len <= headL) return;
    group.add(makeArrow(
      new THREE.Vector3(dir[0], dir[1], dir[2]),
      pos,
      len,
      color,
      headL,
      headW
    ));
  });
  return group;
}

function _rebuildScene() {
  if (!_ellGroup) return;
  while (_ellGroup.children.length) _ellGroup.remove(_ellGroup.children[0]);
  if (!_frames) return;
  const pos = new THREE.Vector3(_frames[6][3], _frames[6][7], _frames[6][11]);
  const r = robotReach * 0.20;
  if (_showTransEl  && _transData) _ellGroup.add(_makeEllipsoid(_transData, pos, 0x00bcd4, r));
  if (_showTransVec && _transData) _ellGroup.add(_makeAxisArrows(_transData, pos, 0x00bcd4, r));
  if (_showRotEl    && _rotData)   _ellGroup.add(_makeEllipsoid(_rotData,   pos, 0xab47bc, r * 0.6));
  if (_showRotVec   && _rotData)   _ellGroup.add(_makeAxisArrows(_rotData,   pos, 0xab47bc, r * 0.6));
}

// Keep old name as alias so any stray internal call still works
const _rebuildEllipsoids = _rebuildScene;

// ─── Public entry points ─────────────────────────────────────────────────────
function computeAndDisplayJacobian(q_deg) {
  if (!currentDH) return;
  _frames   = forwardKin(currentDH, q_deg);
  const q_rad = q_deg.map(v => v * Math.PI / 180);
  _Jspatial = currentJacRows ? evalJacobian(currentJacRows, q_rad) : numericalJacobian(currentDH, q_rad);

  const tv = eig3sym(aat3([_Jspatial[0],_Jspatial[1],_Jspatial[2]]));
  const rv = eig3sym(aat3([_Jspatial[3],_Jspatial[4],_Jspatial[5]]));
  _transData = { sigmas: tv.vals.map(Math.sqrt), vecs: tv.vecs };
  _rotData   = { sigmas: rv.vals.map(Math.sqrt), vecs: rv.vecs };

  _refreshDisplay();
}

function _refreshDisplay() {
  if (!_Jspatial) return;

  _renderJacobianTable(_Jspatial);

  const w = Math.abs(_det6(_Jspatial));
  const [mant, exp] = w.toExponential(3).split('e');
  const formula = katex.renderToString('\\mu = \\sqrt{\\det(J J^{T})} =', {throwOnError: false});
  document.getElementById('manip-overview').innerHTML =
    `<div class="manip-row" style="font-size:0.82rem;align-items:baseline;gap:6px">` +
      `${formula}<span class="m-val">${w.toExponential(3)}</span>` +
    `</div>`;

  _renderManipBlock('manip-trans', 't', _transData);
  _renderManipBlock('manip-rot',   'r', _rotData);
  _syncCheckboxes();
  _rebuildScene();

  document.getElementById('analysis-empty').style.display = 'none';
  document.getElementById('analysis-content').style.display = '';
}
