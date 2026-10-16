'use strict';
/* ======================================================================
 *  Лабораторна робота №1 · Лінійні перетворення · Варіант 4
 *  Корявський Всеволод Андрійович
 *
 *  Архітектура: об'єктно-орієнтована.
 *    Mat3        — матриця 3×3 (вектор-стовпець, p' = M·p), вся арифметика власна
 *    Seg / Arc   — геометричні примітиви деталі (відрізок і дуга кола, крок 1°)
 *    PartDef     — параметрична модель деталі з креслення №4
 *    CoordSystem — звичайна / афінна / проективна система координат
 *    Scene       — рендер у SVG, події миші, перетягування опорної точки
 *  Матриці методички записані рядками (P·T); тут використано транспоновану
 *  форму T' = Tᵀ, результат тотожний.
 * ==================================================================== */

const CELL = 10, GRID_W = 22, GRID_H = 15;
const MARGIN = { x: 64, y: 54 };
const REF_PPC = 35;
const BASE = { x: 20, y: 20 };
const rad = d => d * Math.PI / 180;
const $ = id => document.getElementById(id);

/* ---------------------------------------------------------------- Mat3 */
class Mat3 {
    constructor(a = [1, 0, 0, 0, 1, 0, 0, 0, 1]) { this.a = a; }
    static translate(dx, dy) { return new Mat3([1, 0, dx, 0, 1, dy, 0, 0, 1]); }
    static scale(sx, sy = sx, cx = 0, cy = 0) { return new Mat3([sx, 0, cx * (1 - sx), 0, sy, cy * (1 - sy), 0, 0, 1]); }
    static rotate(deg, cx, cy) {
        const c = Math.cos(rad(deg)), s = Math.sin(rad(deg));
        return new Mat3([c, -s, cx - c * cx + s * cy, s, c, cy - s * cx - c * cy, 0, 0, 1]);
    }
    static pointReflection(cx, cy) { return new Mat3([-1, 0, 2 * cx, 0, -1, 2 * cy, 0, 0, 1]); }
    /** афінна: образи точок (CELL;0), (0;CELL) та початку координат */
    static affine(q) {
        return new Mat3([q.Xx / CELL, q.Yx / CELL, q.Ox, q.Xy / CELL, q.Yy / CELL, q.Oy, 0, 0, 1]);
    }
    /** проективна: кінці осей та ваги (формула (17), транспонована) */
    static projective(q) {
        return new Mat3([q.Xx * q.wX, q.Yx * q.wY, q.Ox * q.w0,
                         q.Xy * q.wX, q.Yy * q.wY, q.Oy * q.w0,
                         q.wX,        q.wY,        q.w0]);
    }
    /** this ∘ other  (спочатку other, потім this) */
    mul(o) {
        const A = this.a, B = o.a, C = new Array(9);
        for (let r = 0; r < 3; r++)
            for (let c = 0; c < 3; c++)
                C[r * 3 + c] = A[r * 3] * B[c] + A[r * 3 + 1] * B[3 + c] + A[r * 3 + 2] * B[6 + c];
        return new Mat3(C);
    }
    inverse() {
        const [a, b, c, d, e, f, g, h, i] = this.a;
        const A = e * i - f * h, B = -(d * i - f * g), C = d * h - e * g;
        const det = a * A + b * B + c * C;
        if (Math.abs(det) < 1e-12) return null;
        return new Mat3([A, -(b * i - c * h), b * f - c * e,
                         B, a * i - c * g, -(a * f - c * d),
                         C, -(a * h - b * g), a * e - b * d].map(v => v / det));
    }
    /** образ точки з діленням на w; null, якщо точка за площиною w ≤ 0 */
    apply(x, y) {
        const m = this.a, w = m[6] * x + m[7] * y + m[8];
        if (w < 1e-6) return null;
        return [(m[0] * x + m[1] * y + m[2]) / w, (m[3] * x + m[4] * y + m[5]) / w];
    }
    /** запис у вигляді рядкової матриці P·T, як у методичці */
    rows() { const m = this.a; return [[m[0], m[3], m[6]], [m[1], m[4], m[7]], [m[2], m[5], m[8]]]; }
}

/* ---------------------------------------------------------- примітиви */
class Seg {
    constructor(a, b, kind = 'edge') { this.a = a; this.b = b; this.kind = kind; }
    sample() { return [this.a, this.b]; }
}
class Arc {
    constructor(cx, cy, r, from, to, kind = 'edge') { Object.assign(this, { cx, cy, r, from, to, kind }); }
    sample() {
        const out = [], steps = Math.max(1, Math.round(Math.abs(this.to - this.from)));
        for (let i = 0; i <= steps; i++) {
            const t = rad(this.from + (this.to - this.from) * i / steps);
            out.push([this.cx + this.r * Math.cos(t), this.cy + this.r * Math.sin(t)]);
        }
        return out;
    }
}
const chain = (names, pts, kind = 'edge') => {
    const out = [];
    for (let i = 0; i + 1 < names.length; i++) out.push(new Seg(pts[names[i]], pts[names[i + 1]], kind));
    return out;
};

/* ------------------------------------------------- деталь (креслення 4) */
class PartDef {
    static FIELDS = [
        ['W', 'Загальна ширина', 100, 160, 125],
        ['H', 'Загальна висота', 60, 110, 80],
        ['gap', 'Просвіт між щоками', 20, 60, 35],
        ['span', 'Відстань між днами вирізів', 60, 110, 80],
        ['open', 'Розкриття V-вирізу', 30, 76, 60],
        ['inner', 'Висота дна вирізу', 6, 40, 20],
        ['pw', 'Ширина пластини', 40, 100, 63],
        ['ph', 'Висота пластини', 40, 76, 60],
        ['dh', 'Діаметр отвору Ø', 10, 55, 30]
    ];
    static defaults() { return Object.fromEntries(PartDef.FIELDS.map(f => [f[0], f[4]])); }

    /** вершини зліва направо: C B A / D E F G H K L / Q M — літери з креслення */
    static build(p) {
        const W = p.W, H = p.H, gap = Math.min(p.gap, p.pw - 4);
        const jaw = (W - gap) / 2;
        const depth = Math.min((W - Math.min(p.span, W - 10)) / 2, jaw - 4);
        const open = Math.min(p.open, H - 4), flat = Math.min(p.inner, open - 2);
        const ph = Math.min(p.ph, H - 2), pw = Math.min(p.pw, W - 4), dh = Math.min(p.dh, ph - 4, pw - 4);
        const mid = H / 2, x0 = (W - pw) / 2, x1 = (W + pw) / 2, yb = mid - ph / 2, yt = mid + ph / 2;
        const V = {
            C: [0, H], B: [jaw, H], A: [jaw, yt], D: [0, mid + open / 2], E: [depth, mid + flat / 2],
            F: [depth, mid - flat / 2], G: [0, mid - open / 2], H: [0, 0], K: [jaw, 0], L: [jaw, yb],
            Q: [x0, yt], M: [x0, yb], N: [x1, yb], P: [x1, yt]
        };
        const mirror = ([x, y]) => [W - x, y];
        for (const k of ['C', 'B', 'A', 'D', 'E', 'F', 'G', 'H', 'K', 'L']) V[k + "'"] = mirror(V[k]);

        const prims = [];
        const half = suffix => chain(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'K', 'L'].map(n => n + suffix), V);
        prims.push(...half(''), ...half("'"));
        prims.push(...chain(['Q', 'P', 'N', 'M', 'Q'], V));
        prims.push(new Arc(W / 2, mid, dh / 2, 0, 360));
        prims.push(new Seg([depth - 6, mid], [W - depth + 6, mid], 'axial'));
        prims.push(new Seg([W / 2, yb - 6], [W / 2, yt + 6], 'axial'));
        return { prims, V, O: [W / 2, mid] };
    }
}

/* -------------------------------------------------- система координат */
class CoordSystem {
    constructor() {
        this.mode = 'std';
        this.aff = { Xx: 20, Xy: 10, Yx: 0, Yy: 10, Ox: 0, Oy: 0 };
        this.prj = { Xx: 520, Xy: 0, wX: 2, Yx: 0, Yy: 430, wY: 3, Ox: 0, Oy: 0, w0: 400 };
    }
    matrix(ppc) {
        const k = ppc / CELL;
        if (this.mode === 'aff') return Mat3.scale(k).mul(Mat3.affine(this.aff));
        if (this.mode === 'prj') return Mat3.scale(ppc / REF_PPC).mul(Mat3.projective(this.prj));
        return Mat3.scale(k);
    }
}

/* ------------------------------------------------------------- стан */
const state = {
    ppc: 35,
    part: PartDef.defaults(),
    pivot: { x: 82.5, y: 60 },
    shift: { x: 0, y: 0 },
    angle: 0,
    mirror: false,
    scale: { x: 1, y: 1 },
    cs: new CoordSystem()
};
const partMatrix = () => {
    const { x, y } = state.pivot;
    let m = Mat3.scale(state.scale.x, state.scale.y, x, y);
    if (state.mirror) m = Mat3.pointReflection(x, y).mul(m);
    m = Mat3.rotate(state.angle, x, y).mul(m);
    return Mat3.translate(state.shift.x, state.shift.y).mul(m);
};

/* ------------------------------------------------------------ Scene */
const SVGNS = 'http://www.w3.org/2000/svg';
class Scene {
    constructor() {
        this.svg = $('scene');
        this.g = { grid: $('gGrid'), part: $('gPart'), marks: $('gMarks'), text: $('gText') };
        this.W = 800; this.H = 600;
        this.resize();
        new ResizeObserver(() => this.resize()).observe($('svgWrap'));
        this.initPointer();
    }
    resize() {
        const r = $('svgWrap').getBoundingClientRect();
        this.W = r.width; this.H = r.height;
        this.svg.setAttribute('viewBox', `0 0 ${this.W} ${this.H}`);
        this.draw();
    }
    /** SVG-координати: Y розвернуто групою (scale(1,-1)), як у Qt */
    flip() { return `matrix(1 0 0 -1 ${MARGIN.x} ${this.H - MARGIN.y})`; }
    toPx(q) { return [q[0] + MARGIN.x, this.H - MARGIN.y - q[1]]; }
    fromPx(sx, sy) { return [sx - MARGIN.x, this.H - MARGIN.y - sy]; }

    el(parent, tag, attrs, text) {
        const e = document.createElementNS(SVGNS, tag);
        for (const k in attrs) e.setAttribute(k, attrs[k]);
        if (text !== undefined) e.textContent = text;
        parent.appendChild(e); return e;
    }
    /** ламана у "d" атрибут, із розривами там, де w ≤ 0 */
    pathData(pts, M, close) {
        let d = '', pen = false;
        const seq = close ? pts.concat([pts[0]]) : pts;
        for (const p of seq) {
            const q = M.apply(p[0], p[1]);
            if (!q) { pen = false; continue; }
            d += (pen ? 'L' : 'M') + q[0].toFixed(2) + ' ' + q[1].toFixed(2);
            pen = true;
        }
        return d;
    }
    lerp(a, b, n = 36) { return Array.from({ length: n + 1 }, (_, i) => [a[0] + (b[0] - a[0]) * i / n, a[1] + (b[1] - a[1]) * i / n]); }

    draw() {
        const Tc = state.cs.matrix(state.ppc), Tm = partMatrix();
        for (const k in this.g) this.g[k].replaceChildren();
        for (const k of ['grid', 'part']) this.g[k].setAttribute('transform', this.flip());

        // сітка та осі (єдиний об'єкт "осі+сітка")
        for (let i = 0; i <= GRID_W; i++) this.el(this.g.grid, 'path', { class: 'gridline', d: this.pathData(this.lerp([i * CELL, 0], [i * CELL, GRID_H * CELL]), Tc) });
        for (let j = 0; j <= GRID_H; j++) this.el(this.g.grid, 'path', { class: 'gridline', d: this.pathData(this.lerp([0, j * CELL], [GRID_W * CELL, j * CELL]), Tc) });
        this.el(this.g.grid, 'path', { class: 'axisline', d: this.pathData(this.lerp([0, 0], [GRID_W * CELL, 0]), Tc) });
        this.el(this.g.grid, 'path', { class: 'axisline', d: this.pathData(this.lerp([0, 0], [0, GRID_H * CELL]), Tc) });
        this.axisDecor(Tc);

        // деталь: T_деталі · T_координат
        const { prims, V, O } = PartDef.build(state.part);
        const full = Tc.mul(Tm).mul(Mat3.translate(BASE.x, BASE.y));   // деталь стоїть у (BASE) на сітці
        const axialOn = $('optAxial').checked;
        for (const pr of prims) {
            if (pr.kind === 'axial' && !axialOn) continue;
            const closed = pr instanceof Arc && Math.abs(pr.to - pr.from) >= 360;
            this.el(this.g.part, 'path', { class: pr.kind, d: this.pathData(pr.sample(), full, closed) });
        }
        if ($('optLabels').checked) this.vertexLabels(V, O, full);

        // опорна точка (обертання/симетрія/масштаб)
        const pv = Tc.apply(state.pivot.x, state.pivot.y);
        if (pv) {
            const [x, y] = this.toPx(pv);
            const grp = this.el(this.g.marks, 'g', { class: 'pivot', id: 'pivot' });
            this.el(grp, 'circle', { cx: x, cy: y, r: 13, fill: 'rgba(209,67,91,.10)', stroke: '#d1435b', 'stroke-width': 1.4 });
            this.el(grp, 'path', { d: `M${x - 9} ${y}H${x + 9}M${x} ${y - 9}V${y + 9}`, stroke: '#d1435b', 'stroke-width': 1.4 });
        }
        $('mPart').textContent = this.fmt(Tm.rows());
        $('mCs').textContent = this.fmt(Tc.rows());
    }
    fmt(rows) { return rows.map(r => r.map(v => (Math.round(v * 1000) / 1000).toFixed(3).padStart(9)).join('')).join('\n'); }

    axisDecor(Tc) {
        const T = this.g.text;
        const tip = (p, q) => {
            const a = Tc.apply(...p), b = Tc.apply(...q); if (!a || !b) return null;
            const A = this.toPx(a), B = this.toPx(b), L = Math.hypot(A[0] - B[0], A[1] - B[1]) || 1;
            return { A, u: [(A[0] - B[0]) / L, (A[1] - B[1]) / L] };
        };
        for (const [name, p, q] of [['X', [GRID_W * CELL, 0], [(GRID_W - .5) * CELL, 0]], ['Y', [0, GRID_H * CELL], [0, (GRID_H - .5) * CELL]]]) {
            const t = tip(p, q); if (!t) continue;
            const { A, u } = t;
            this.el(T, 'path', { d: `M${A[0]} ${A[1]}l${-12 * u[0] + 5 * u[1]} ${-12 * u[1] - 5 * u[0]}M${A[0]} ${A[1]}l${-12 * u[0] - 5 * u[1]} ${-12 * u[1] + 5 * u[0]}`, stroke: '#22252e', 'stroke-width': 1.5, fill: 'none' });
            this.el(T, 'text', { class: 'lbl', x: A[0] + 8, y: A[1] + (name === 'X' ? 5 : -8) }, name.toLowerCase());
        }
        let prev = -99;
        for (let i = 0; i <= GRID_W; i++) {
            const q = Tc.apply(i * CELL, 0); if (!q) continue;
            const s = this.toPx(q); if (s[0] - prev < 24) continue; prev = s[0];
            this.el(T, 'text', { class: 'num', x: s[0] - 6, y: s[1] + 14 }, String(i * CELL));
        }
        prev = 1e9;
        for (let j = 1; j <= GRID_H; j++) {
            const q = Tc.apply(0, j * CELL); if (!q) continue;
            const s = this.toPx(q); if (prev - s[1] < 14) continue; prev = s[1];
            this.el(T, 'text', { class: 'num', x: s[0] - 24, y: s[1] + 3 }, String(j * CELL));
        }
    }
    vertexLabels(V, O, M) {
        const place = (name, p, dx, dy) => {
            const q = M.apply(p[0], p[1]); if (!q) return;
            const s = this.toPx(q);
            this.el(this.g.text, 'circle', { cx: s[0], cy: s[1], r: 2.6, fill: '#22252e' });
            this.el(this.g.text, 'text', { class: 'lbl', x: s[0] + dx, y: s[1] + dy }, name);
        };
        for (const k in V) {
            const top = ['C', 'B', 'A', 'Q', 'P'].includes(k[0]);
            place(k, V[k], k[0] === 'E' || k[0] === 'F' ? 8 : -14, top ? -8 : 16);
        }
        place('O', O, 6, 15);
    }

    initPointer() {
        const svg = this.svg;
        const toWorld = e => {
            const r = svg.getBoundingClientRect();
            const inv = state.cs.matrix(state.ppc).inverse();
            const p = this.fromPx(e.clientX - r.left, e.clientY - r.top);
            return inv ? inv.apply(p[0], p[1]) : null;
        };
        let drag = false;
        svg.addEventListener('pointerdown', e => {
            if (e.target.closest('#pivot')) { drag = true; svg.setPointerCapture(e.pointerId); }
        });
        svg.addEventListener('pointerup', () => { drag = false; });
        svg.addEventListener('pointermove', e => {
            const w = toWorld(e);
            $('hud').textContent = w ? `курсор: (${w[0].toFixed(1)}; ${w[1].toFixed(1)}) мм` : 'курсор: –';
            if (drag && w) {
                state.pivot.x = Math.round(w[0] * 2) / 2; state.pivot.y = Math.round(w[1] * 2) / 2;
                widgets.pvx.set(state.pivot.x); widgets.pvy.set(state.pivot.y);
                this.draw();
            }
        });
        svg.addEventListener('wheel', e => {
            e.preventDefault();
            state.ppc = Math.min(80, Math.max(12, state.ppc + (e.deltaY < 0 ? 2 : -2)));
            widgets.ppc.set(state.ppc); this.draw();
        }, { passive: false });
    }
}

/* ------------------------------------------------------- віджети */
const widgets = {};
/** поле «назва + число + повзунок» */
function field(host, label, min, max, step, value, onChange) {
    const wrap = document.createElement('div'); wrap.className = 'fld';
    wrap.innerHTML = `<span class="n">${label}</span><input type="number" min="${min}" max="${max}" step="${step}"><input type="range" min="${min}" max="${max}" step="${step}">`;
    const [num, rng] = wrap.querySelectorAll('input');
    const set = v => { num.value = v; rng.value = v; };
    const fire = v => { if (Number.isFinite(v)) { set(v); onChange(v); } };
    num.addEventListener('input', () => fire(parseFloat(num.value)));
    rng.addEventListener('input', () => fire(parseFloat(rng.value)));
    set(value); host.appendChild(wrap);
    return { set, el: wrap };
}
function plainField(host, label, key, obj, value) {
    const l = document.createElement('label');
    l.innerHTML = `${label}<input type="number" step="any" value="${value}">`;
    l.querySelector('input').addEventListener('input', e => { obj[key] = parseFloat(e.target.value) || 0; });
    host.appendChild(l);
    return l.querySelector('input');
}

let scene;
function buildUI() {
    widgets.ppc = field($('f_ppc'), 'Пікселів на сантиметр', 12, 80, 1, state.ppc, v => { state.ppc = v; scene.draw(); });

    const part = $('f_part'); part.classList.add('body', 'duo');
    for (const [key, label, min, max, def] of PartDef.FIELDS)
        widgets['p_' + key] = field(part, label, min, max, 1, def, v => { state.part[key] = v; scene.draw(); });

    const duo = id => { const h = $(id); h.className = 'duo'; return h; };
    widgets.pvx = field(duo('f_pivot'), 'x, мм', -50, 220, .5, state.pivot.x, v => { state.pivot.x = v; scene.draw(); });
    widgets.pvy = field($('f_pivot'), 'y, мм', -50, 150, .5, state.pivot.y, v => { state.pivot.y = v; scene.draw(); });
    widgets.shx = field(duo('f_shift'), 'dx', -100, 100, 1, 0, v => { state.shift.x = v; scene.draw(); });
    widgets.shy = field($('f_shift'), 'dy', -100, 100, 1, 0, v => { state.shift.y = v; scene.draw(); });
    widgets.rot = field($('f_rot'), 'Кут, °', -360, 360, 1, 0, v => { state.angle = v; scene.draw(); });
    widgets.scx = field(duo('f_scale'), 'Sx', -3, 3, .05, 1, v => { state.scale.x = v; scene.draw(); });
    widgets.scy = field($('f_scale'), 'Sy', -3, 3, .05, 1, v => { state.scale.y = v; scene.draw(); });

    $('optMirror').onchange = e => { state.mirror = e.target.checked; scene.draw(); };
    $('optLabels').onchange = $('optAxial').onchange = () => scene.draw();

    // афінна та проективна форми
    const cs = state.cs;
    const affIn = {}, prjIn = {};
    for (const [k, t] of [['Xx', 'Xx'], ['Xy', 'Xy'], ['Yx', 'Yx'], ['Yy', 'Yy'], ['Ox', '0x'], ['Oy', '0y']]) affIn[k] = plainField($('f_aff'), t, k, {}, cs.aff[k]);
    for (const [k, t] of [['Xx', 'Xx'], ['Xy', 'Xy'], ['wX', 'wX'], ['Yx', 'Yx'], ['Yy', 'Yy'], ['wY', 'wY'], ['Ox', '0x'], ['Oy', '0y'], ['w0', 'w0']]) prjIn[k] = plainField($('f_prj'), t, k, {}, cs.prj[k]);
    const read = box => Object.fromEntries(Object.entries(box).map(([k, el]) => [k, parseFloat(el.value) || 0]));
    const setMode = mode => {
        cs.mode = mode;
        document.querySelectorAll('#csSeg button').forEach(b => b.classList.toggle('on', b.dataset.cs === mode));
        $('boxAff').hidden = mode !== 'aff'; $('boxPrj').hidden = mode !== 'prj';
        scene.draw();
    };
    document.querySelectorAll('#csSeg button').forEach(b => b.onclick = () => setMode(b.dataset.cs));
    document.querySelectorAll('[data-apply]').forEach(b => b.onclick = () => {
        if (b.dataset.apply === 'aff') cs.aff = read(affIn); else cs.prj = read(prjIn);
        setMode(b.dataset.apply);
    });
    window.setMode = setMode;

    $('btnDefault').onclick = () => {
        state.part = PartDef.defaults();
        for (const [key, , , , def] of PartDef.FIELDS) widgets['p_' + key].set(def);
        scene.draw();
    };
    $('btnReset').onclick = () => {
        Object.assign(state, { shift: { x: 0, y: 0 }, angle: 0, mirror: false, scale: { x: 1, y: 1 } });
        for (const k of ['shx', 'shy', 'rot']) widgets[k].set(0);
        widgets.scx.set(1); widgets.scy.set(1); $('optMirror').checked = false;
        setMode('std');
    };
}

/* початковий стан з адресного рядка (?s={json}) — для знімків екрана */
function applyQuery() {
    const q = new URLSearchParams(location.search).get('s'); if (!q) return;
    try {
        const s = JSON.parse(q);
        if (s.angle !== undefined) { state.angle = s.angle; widgets.rot.set(s.angle); }
        if (s.shift) { Object.assign(state.shift, s.shift); widgets.shx.set(state.shift.x); widgets.shy.set(state.shift.y); }
        if (s.mirror) { state.mirror = true; $('optMirror').checked = true; }
        if (s.scale) { Object.assign(state.scale, s.scale); widgets.scx.set(state.scale.x); widgets.scy.set(state.scale.y); }
        if (s.part) { Object.assign(state.part, s.part); for (const k in s.part) widgets['p_' + k].set(s.part[k]); }
        if (s.pivot) { Object.assign(state.pivot, s.pivot); widgets.pvx.set(state.pivot.x); widgets.pvy.set(state.pivot.y); }
        if (s.mode) window.setMode(s.mode);
    } catch (e) { /* ignore */ }
}

scene = new Scene();
buildUI();
applyQuery();
scene.draw();
