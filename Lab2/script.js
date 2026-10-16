'use strict';
/* ======================================================================
 *  Лабораторна робота №2 · Побудова чудових кривих · Варіант 4
 *  Трисектриса Лоншама   r = a / cos 3φ        Корявський В. А.
 *
 *  Крива задана у ПОЛЯРНИХ координатах (формули (5)-(6) методички), тож
 *  усі характеристики виводяться через r(φ), r′(φ), r″(φ):
 *      довжина дуги   L = ∫ √(r² + r′²) dφ                    (20)
 *      площа сектора  S = ½ ∫ r² dφ
 *      кривина        k = (r² + 2r′² − r·r″) / (r² + r′²)^{3/2}
 *  Об'єктна архітектура: Mat3, Quadrature, Roots, Trisectrix, Scene, App.
 *  SVG малює у СВІТОВИХ одиницях (група з матрицею масштабу), а над
 *  кривою всі евклідові перетворення виконуються власним класом Mat3.
 * ==================================================================== */

const PI = Math.PI, SQ3 = Math.sqrt(3);
const rad = d => d * PI / 180, deg = r => r * 180 / PI;
const $ = id => document.getElementById(id);

/* ------------------------------------------------------------ Mat3 */
class Mat3 {
    constructor(m = [1, 0, 0, 0, 1, 0, 0, 0, 1]) { this.m = m; }
    static translate(dx, dy) { return new Mat3([1, 0, dx, 0, 1, dy, 0, 0, 1]); }
    static rotate(d, cx, cy) {
        const c = Math.cos(rad(d)), s = Math.sin(rad(d));
        return new Mat3([c, -s, cx - c * cx + s * cy, s, c, cy - s * cx - c * cy, 0, 0, 1]);
    }
    /** this ∘ o (спочатку o) */
    mul(o) {
        const a = this.m, b = o.m, r = new Array(9);
        for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) r[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
        return new Mat3(r);
    }
    apply([x, y]) { const m = this.m; return [m[0] * x + m[1] * y + m[2], m[3] * x + m[4] * y + m[5]]; }
    applyDir([x, y]) { const m = this.m; return [m[0] * x + m[1] * y, m[3] * x + m[4] * y]; }
}

/* --------------------------------------------- чисельні методи */
class Quadrature {
    static X = [0, 0.5384693101056831, -0.5384693101056831, 0.9061798459386640, -0.9061798459386640];
    static Wt = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
    /** складена 5-вузлова квадратура Гауса–Лежандра */
    static gauss(f, lo, hi, panels = 160) {
        const h = (hi - lo) / panels; let sum = 0;
        for (let p = 0; p < panels; p++) {
            const mid = lo + (p + 0.5) * h;
            for (let k = 0; k < 5; k++) sum += Quadrature.Wt[k] * f(mid + Quadrature.X[k] * h / 2);
        }
        return sum * h / 2;
    }
}
class Roots {
    static bisect(f, lo, hi, it = 70) {
        let flo = f(lo);
        for (let i = 0; i < it; i++) {
            const mid = (lo + hi) / 2, fm = f(mid);
            if (flo * fm <= 0) hi = mid; else { lo = mid; flo = fm; }
        }
        return (lo + hi) / 2;
    }
    static scan(f, lo, hi, n = 2500) {
        const out = []; let x0 = lo, f0 = f(x0);
        for (let i = 1; i <= n; i++) {
            const x1 = lo + (hi - lo) * i / n, f1 = f(x1);
            if (Number.isFinite(f0) && Number.isFinite(f1) && f0 * f1 < 0) out.push(Roots.bisect(f, x0, x1));
            x0 = x1; f0 = f1;
        }
        return out;
    }
}

/* ------------------------------------------------------------ крива */
class Trisectrix {
    static BRANCHES = [[-PI / 6, PI / 6], [PI / 6, PI / 2], [PI / 2, 5 * PI / 6]];
    static POLES = [PI / 6, PI / 2, 5 * PI / 6];
    constructor(a) { this.a = a; }

    r(p) { return this.a / Math.cos(3 * p); }
    dr(p) { const c = Math.cos(3 * p); return 3 * this.a * Math.sin(3 * p) / (c * c); }
    d2r(p) { const c = Math.cos(3 * p), s = Math.sin(3 * p); return 9 * this.a * (1 + s * s) / (c * c * c); }

    xy(p) { const r = this.r(p); return [r * Math.cos(p), r * Math.sin(p)]; }
    d1(p) { const r = this.r(p), r1 = this.dr(p), c = Math.cos(p), s = Math.sin(p); return [r1 * c - r * s, r1 * s + r * c]; }
    arcElement(p) { return Math.hypot(this.r(p), this.dr(p)); }
    kNumerator(p) { const r = this.r(p); return r * r + 2 * this.dr(p) ** 2 - r * this.d2r(p); }
    curvature(p) { return this.kNumerator(p) / this.arcElement(p) ** 3; }

    hasPole(p1, p2) { return Trisectrix.POLES.some(q => q > Math.min(p1, p2) && q < Math.max(p1, p2)); }
    arcLength(p1, p2) { return this.hasPole(p1, p2) ? Infinity : Quadrature.gauss(p => this.arcElement(p), p1, p2); }
    sectorArea(p1, p2) { return this.hasPole(p1, p2) ? Infinity : 0.5 * Quadrature.gauss(p => this.r(p) ** 2, p1, p2); }
    inflections() {
        return Trisectrix.BRANCHES.flatMap(([lo, hi]) => Roots.scan(p => this.kNumerator(p), lo + 2e-3, hi - 2e-3));
    }
    /** три асимптоти n·P = −a/3, n = (cosθ, sinθ), θ = 0°, 120°, 240° */
    asymptotes() { return [0, 120, 240].map(th => ({ th, n: [Math.cos(rad(th)), Math.sin(rad(th))], c: -this.a / 3 })); }
    triangleArea() { return SQ3 * this.a * this.a / 3; }
}

/* ------------------------------------------------------------ стан */
const state = {
    a: 2, sc: 40, phi: 12,
    shift: [0, 0], rot: { x: 0, y: 0, a: 0 },
    view: { ox: 380, oy: 300 },
    show: { tan: true, nor: true, asy: true, circ: false },
    calc: [-15, 25],
    samples: [], anim: { on: false, t0: 0 }
};
const eucl = () => Mat3.translate(state.shift[0], state.shift[1]).mul(Mat3.rotate(state.rot.a, state.rot.x, state.rot.y));

/* ------------------------------------------------------------ Scene */
const NS = 'http://www.w3.org/2000/svg';
class Scene {
    constructor() {
        this.svg = $('scene');
        this.W = 900; this.H = 560;
        this.world = this.mk('g', {}, this.svg);
        this.ui = this.mk('g', {}, this.svg);
        this.gridG = this.mk('g', {}, this.world);
        this.curveG = this.mk('g', {}, this.world);
        this.itemG = this.mk('g', {}, this.world);
        this.firstFit = true;
        new ResizeObserver(() => this.fit()).observe($('svgWrap'));
        this.pointer();
    }
    mk(tag, attrs = {}, parent) {
        const e = document.createElementNS(NS, tag);
        for (const k in attrs) e.setAttribute(k, attrs[k]);
        if (parent) parent.appendChild(e);
        return e;
    }
    fit() {
        const r = $('svgWrap').getBoundingClientRect();
        this.W = r.width; this.H = r.height;
        this.svg.setAttribute('viewBox', `0 0 ${this.W} ${this.H}`);
        if (this.firstFit) { state.view.ox = this.W * 0.38; state.view.oy = this.H * 0.5; this.firstFit = false; }
        app.draw();
    }
    toScreen([x, y]) { return [state.view.ox + x * state.sc, state.view.oy - y * state.sc]; }
    toWorld(sx, sy) { return [(sx - state.view.ox) / state.sc, (state.view.oy - sy) / state.sc]; }

    path(pts, cls, color) {
        let d = '', pen = false;
        for (const p of pts) {
            if (Math.abs(p[0]) > 600 || Math.abs(p[1]) > 600 || !Number.isFinite(p[0] + p[1])) { pen = false; continue; }
            d += (pen ? 'L' : 'M') + p[0].toFixed(4) + ' ' + p[1].toFixed(4); pen = true;
        }
        const el = this.mk('path', { d, class: cls }, this.curveG);
        if (color) el.setAttribute('stroke', color);
        return el;
    }
    line(p, q, cls, parent = this.itemG) { return this.mk('line', { x1: p[0], y1: p[1], x2: q[0], y2: q[1], class: cls }, parent); }

    drawGrid() {
        this.gridG.replaceChildren(); this.ui.replaceChildren();
        const [x0, y1] = this.toWorld(0, 0), [x1, y0] = this.toWorld(this.W, this.H);
        for (let i = Math.floor(x0); i <= Math.ceil(x1); i++) this.line([i, y0], [i, y1], i % 5 ? 'gridl' : 'gridM', this.gridG);
        for (let j = Math.floor(y0); j <= Math.ceil(y1); j++) this.line([x0, j], [x1, j], j % 5 ? 'gridl' : 'gridM', this.gridG);
        this.line([x0, 0], [x1, 0], 'axisX', this.gridG); this.line([0, y0], [0, y1], 'axisY', this.gridG);
        const o = this.toScreen([0, 0]);
        for (let i = Math.ceil(x0); i <= x1; i++) if (i && i % 2 === 0) this.text(i, o[0] + i * state.sc - 4, o[1] + 15, 'lbl');
        for (let j = Math.ceil(y0); j <= y1; j++) if (j && j % 2 === 0) this.text(j, o[0] + 7, o[1] - j * state.sc + 4, 'lbl');
        this.text('x', this.W - 16, o[1] - 8, 'mk'); this.text('y', o[0] + 9, 16, 'mk'); this.text('O', o[0] - 14, o[1] + 15, 'mk');
    }
    text(s, x, y, cls) { this.mk('text', { x, y, class: cls }, this.ui).textContent = s; }

    pointer() {
        let drag = null;
        this.svg.addEventListener('pointerdown', e => {
            drag = { x: e.clientX, y: e.clientY, ox: state.view.ox, oy: state.view.oy, moved: 0 };
            this.svg.setPointerCapture(e.pointerId); this.svg.classList.add('drag');
        });
        this.svg.addEventListener('pointermove', e => {
            if (!drag) return;
            const dx = e.clientX - drag.x, dy = e.clientY - drag.y; drag.moved = Math.max(drag.moved, Math.hypot(dx, dy));
            state.view.ox = drag.ox + dx; state.view.oy = drag.oy + dy; app.draw();
        });
        this.svg.addEventListener('pointerup', e => {
            this.svg.classList.remove('drag');
            if (drag && drag.moved < 4) {
                const r = this.svg.getBoundingClientRect();
                app.pickAt(e.clientX - r.left, e.clientY - r.top);
            }
            drag = null;
        });
        this.svg.addEventListener('wheel', e => {
            e.preventDefault();
            state.sc = Math.min(90, Math.max(10, state.sc + (e.deltaY < 0 ? 2 : -2)));
            app.syncInputs(); app.draw();
        }, { passive: false });
    }
}

/* -------------------------------------------------------------- App */
const f4 = v => Number.isFinite(v) ? (Math.round(v * 1e4) / 1e4).toString() : '∞';
const f5 = v => Number.isFinite(v) ? v.toFixed(5) : '∞';
const BRANCH_COLORS = ['#22252e', '#4f46e5', '#c026d3'];

const app = {
    curve: new Trisectrix(2),
    draw() {
        const c = this.curve = new Trisectrix(state.a), T = eucl(), sc = scene;
        $('scene').style.setProperty('--sc', state.sc);
        sc.world.setAttribute('transform', `matrix(${state.sc} 0 0 ${-state.sc} ${state.view.ox} ${state.view.oy})`);
        sc.drawGrid();
        sc.curveG.replaceChildren(); sc.itemG.replaceChildren();

        // побудова кривої: цикл за φ, зберігаємо у контейнері кожну 20-ту точку
        state.samples = [];
        const N = 1400;
        Trisectrix.BRANCHES.forEach(([lo, hi], bi) => {
            const eps = 1e-4, pts = [];
            for (let i = 0; i <= N; i++) {
                const u = i / N, phi = lo + eps + (hi - lo - 2 * eps) * (1 - Math.cos(PI * u)) / 2;
                const w = c.xy(phi), p = T.apply(w);
                pts.push(p);
                if (i % 20 === 0) state.samples.push({ phi, w, p, br: bi });
            }
            sc.path(pts, 'br', BRANCH_COLORS[bi]);
        });
        state.samples.sort((u, v) => u.phi - v.phi);

        const L = 700;
        if (state.show.asy) {
            const tri = c.asymptotes().map(l => T.apply([-l.n[1] * 0 + 2 * c.a / 3 * l.n[0], 2 * c.a / 3 * l.n[1]]));
            sc.mk('polygon', { points: tri.map(p => p.join(',')).join(' '), class: 'tri' }, sc.itemG);
            for (const l of c.asymptotes()) {
                const p0 = [l.n[0] * l.c, l.n[1] * l.c], d = [-l.n[1], l.n[0]];
                sc.line(T.apply([p0[0] - d[0] * L, p0[1] - d[1] * L]), T.apply([p0[0] + d[0] * L, p0[1] + d[1] * L]), 'asy');
            }
        }
        this.drawPoint(c, T);
        this.updateResults(c, T);
        if (!$('dTab').hidden) this.fillTable(T);
    },
    pointInfo(c, T) {
        const phi = rad(state.phi), w = c.xy(phi), d = c.d1(phi), v = Math.hypot(d[0], d[1]);
        const k = c.curvature(phi);
        return { phi, w, p: T.apply(w), tau: [d[0] / v, d[1] / v], k, R: 1 / Math.abs(k), r: c.r(phi) };
    },
    drawPoint(c, T) {
        const m = this.pointInfo(c, T), half = 230 / state.sc, g = scene.itemG;
        const seg = (dir, cls) => scene.line(T.apply([m.w[0] - dir[0] * half, m.w[1] - dir[1] * half]), T.apply([m.w[0] + dir[0] * half, m.w[1] + dir[1] * half]), cls);
        if (state.show.tan) seg(m.tau, 'tan');
        if (state.show.nor) seg([-m.tau[1], m.tau[0]], 'nor');
        if (state.show.circ && m.R < 300) {
            const s = Math.sign(m.k) || 1, ctr = T.apply([m.w[0] - s * m.tau[1] * m.R, m.w[1] + s * m.tau[0] * m.R]);
            scene.mk('circle', { cx: ctr[0], cy: ctr[1], r: m.R, class: 'circ' }, g);
        }
        // маркер M — у екранних координатах, щоб радіус не залежав від масштабу
        const s = scene.toScreen(m.p);
        scene.mk('circle', { cx: s[0], cy: s[1], r: 6, fill: '#f59e0b', stroke: '#fff', 'stroke-width': 2 }, scene.ui);
        scene.text('M', s[0] + 10, s[1] - 9, 'mk');
        if (state.rot.a) { const q = scene.toScreen([state.rot.x, state.rot.y]); scene.mk('path', { d: `M${q[0] - 9} ${q[1]}H${q[0] + 9}M${q[0]} ${q[1] - 9}V${q[1] + 9}`, stroke: '#d1435b', 'stroke-width': 1.4 }, scene.ui); }
    },

    /* ----- дані нижньої панелі ----- */
    buildResultCards() {
        $('dRes').innerHTML = `
        <div class="stat"><small>Довжина дуги L</small><b id="vL"></b>
            <div class="row">φ₁<input id="c1" type="number" value="-15" step="1">φ₂<input id="c2" type="number" value="25" step="1"></div></div>
        <div class="stat"><small>Площа сектора ½∫r²dφ</small><b id="vS"></b><p>на тому ж проміжку</p></div>
        <div class="stat"><small>Площа трикутника з асимптот</small><b id="vT"></b><p>√3·a²/3; між кривою й асимптотою площа розбіжна (∞)</p></div>
        <div class="stat"><small>Точка M</small><b id="vM"></b><p id="vM2"></p></div>
        <div class="stat"><small>Радіус кривини R у точці M</small><b id="vR"></b><p id="vR2"></p></div>
        <div class="stat"><small>Точки перегину</small><b id="vI"></b><p id="vI2"></p></div>
        <div class="stat"><small>Асимптоти</small><p id="vA" style="margin-top:0;color:var(--ink)"></p></div>
        <div class="stat"><small>Рівняння дотичної / нормалі</small><p id="vE" style="margin-top:0;color:var(--ink)"></p></div>`;
        $('c1').oninput = $('c2').oninput = () => { state.calc = [parseFloat($('c1').value) || 0, parseFloat($('c2').value) || 0]; this.draw(); };
    },
    updateResults(c, T) {
        const m = this.pointInfo(c, T), [p1, p2] = state.calc.map(rad);
        const L = c.arcLength(p1, p2), S = c.sectorArea(p1, p2);
        $('vL').textContent = Number.isFinite(L) ? f5(L) : '∞ (полюс)';
        $('vS').textContent = Number.isFinite(S) ? f5(S) : '∞ (полюс)';
        $('vT').textContent = f5(c.triangleArea());
        $('vM').textContent = `(${f4(m.p[0])}; ${f4(m.p[1])})`;
        $('vM2').textContent = `φ = ${f4(state.phi)}°, r = ${f4(m.r)}; до T: (${f4(m.w[0])}; ${f4(m.w[1])})`;
        $('vR').textContent = f5(m.R);
        $('vR2').textContent = `k = ${f5(m.k)};  у вершині (φ=0): R = a/8 = ${f5(c.a / 8)}`;
        const inf = c.inflections();
        $('vI').textContent = inf.length ? inf.map(p => f4(deg(p)) + '°').join(', ') : 'немає';
        $('vI2').textContent = inf.length ? '' : 'r² + 2r′² − r·r″ = −8a²/cos²3φ ≠ 0';
        $('vA').innerHTML = c.asymptotes().map(l => `${f4(l.n[0])}·x + ${f4(l.n[1])}·y = ${f4(l.c)}`).join('<br>');
        const t = m.tau, d = Math.abs(t[0]) < 1e-9 ? null : t[1] / t[0];
        $('vE').innerHTML = d === null ? `x = ${f4(m.w[0])}<br>y = ${f4(m.w[1])}`
            : `y = ${f4(m.w[1])} ${d < 0 ? '−' : '+'} ${f4(Math.abs(d))}·(x ${m.w[0] < 0 ? '+' : '−'} ${f4(Math.abs(m.w[0]))})<br>`
              + `y = ${f4(m.w[1])} ${(-1 / d) < 0 ? '−' : '+'} ${f4(Math.abs(1 / d))}·(x ${m.w[0] < 0 ? '+' : '−'} ${f4(Math.abs(m.w[0]))})`;
    },
    fillTable(T) {
        const rows = state.samples.map((s, i) => `<tr data-i="${i}"><td>${i + 1}</td><td>${deg(s.phi).toFixed(2)}</td><td>${s.w[0].toFixed(3)}</td><td>${s.w[1].toFixed(3)}</td><td>${s.p[0].toFixed(3)}</td><td>${s.p[1].toFixed(3)}</td></tr>`);
        $('tbl').tBodies[0].innerHTML = rows.join('');
        const near = this.nearestSample(); const tr = $('tbl').tBodies[0].rows[near]; if (tr) tr.classList.add('sel');
    },
    nearestSample() {
        let best = 0, bd = 1e9;
        state.samples.forEach((s, i) => { const d = Math.abs(deg(s.phi) - state.phi); if (d < bd) { bd = d; best = i; } });
        return best;
    },
    step(dir) {
        const i = Math.min(state.samples.length - 1, Math.max(0, this.nearestSample() + dir));
        this.setPhi(deg(state.samples[i].phi));
    },
    pickAt(sx, sy) {
        let best = null, bd = 30;
        for (const s of state.samples) {
            const q = scene.toScreen(s.p), d = Math.hypot(q[0] - sx, q[1] - sy);
            if (d < bd) { bd = d; best = s; }
        }
        if (best) this.setPhi(deg(best.phi));
    },
    setPhi(p) { state.phi = Math.round(p * 100) / 100; this.syncInputs(); this.draw(); },

    syncInputs() {
        $('a').value = $('aR').value = state.a; $('sc').value = $('scR').value = state.sc;
        $('phi').value = $('phiR').value = state.phi;
        $('shx').value = state.shift[0]; $('shy').value = state.shift[1];
        $('rx').value = state.rot.x; $('ry').value = state.rot.y; $('ra').value = $('raR').value = state.rot.a;
    },

    /* ----- анімація ----- */
    tick(now) {
        if (!state.anim.on) return;
        const t = (now - state.anim.t0) / 1000, v = +$('anV').value * 0.25;
        if ($('anA').checked) state.a = Math.round((2.5 + 1.5 * Math.sin(t * v)) * 100) / 100;
        if ($('anM').checked) state.phi = Math.round(25 * Math.sin(t * v * 1.7) * 100) / 100;
        if ($('anR').checked) state.rot.a = Math.round(((t * v * 25) % 360) * 10) / 10;
        this.syncInputs(); this.draw();
        requestAnimationFrame(n => this.tick(n));
    },
    bind() {
        const num = id => { const v = parseFloat($(id).value); return Number.isFinite(v) ? v : 0; };
        const link = (id, idR, apply) => { for (const k of [id, idR]) $(k).oninput = () => { apply(num(k)); this.syncInputs(); this.draw(); }; };
        link('a', 'aR', v => state.a = Math.min(6, Math.max(0.5, v)));
        link('sc', 'scR', v => state.sc = Math.min(90, Math.max(10, v)));
        link('phi', 'phiR', v => state.phi = v);
        link('ra', 'raR', v => state.rot.a = v);
        for (const id of ['shx', 'shy', 'rx', 'ry']) $(id).oninput = () => {
            state.shift = [num('shx'), num('shy')]; state.rot.x = num('rx'); state.rot.y = num('ry'); this.draw();
        };
        $('btnDefault').onclick = () => { state.a = 2; state.sc = 40; state.phi = 12; this.syncInputs(); this.draw(); };
        $('btnReset').onclick = () => { state.shift = [0, 0]; state.rot = { x: 0, y: 0, a: 0 }; this.syncInputs(); this.draw(); };
        for (const [id, key] of [['vTan', 'tan'], ['vNor', 'nor'], ['vAsy', 'asy'], ['vCirc', 'circ']]) $(id).onchange = e => { state.show[key] = e.target.checked; this.draw(); };
        $('pPrev').onclick = () => this.step(-1); $('pNext').onclick = () => this.step(1);
        $('btnAnim').onclick = () => {
            state.anim.on = !state.anim.on; $('btnAnim').textContent = state.anim.on ? '■ Зупинити' : '▶ Запустити';
            if (state.anim.on) { state.anim.t0 = performance.now(); requestAnimationFrame(n => this.tick(n)); }
        };
        document.querySelectorAll('#dockTabs button').forEach(b => b.onclick = () => {
            document.querySelectorAll('#dockTabs button').forEach(x => x.classList.toggle('on', x === b));
            $('dRes').hidden = b.dataset.d !== 'res'; $('dTab').hidden = b.dataset.d !== 'tab';
            this.draw();
        });
        $('tbl').tBodies[0].onclick = e => { const tr = e.target.closest('tr'); if (tr) this.setPhi(deg(state.samples[+tr.dataset.i].phi)); };
    },
    applyQuery() {
        const q = new URLSearchParams(location.search).get('s'); if (!q) return;
        try {
            const s = JSON.parse(q);
            for (const k of ['a', 'phi', 'sc']) if (s[k] !== undefined) state[k] = s[k];
            if (s.shift) state.shift = s.shift;
            if (s.rot) Object.assign(state.rot, s.rot);
            if (s.circ) { state.show.circ = true; $('vCirc').checked = true; }
            if (s.tab) document.querySelector(`#dockTabs button[data-d="${s.tab}"]`).click();
            if (s.anim) $('btnAnim').click();
        } catch (e) { /* ignore */ }
    }
};

const scene = new Scene();
app.buildResultCards(); app.bind(); app.applyQuery(); app.syncInputs(); scene.fit();
