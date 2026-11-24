'use strict';
/* ======================================================================
 *  Лабораторна робота №3 · Моделювання криволінійних контурів · Варіант 4
 *  Контур №27 «Скакалка», метод №3 — криві другого порядку, задані
 *  в ІНЖЕНЕРНОМУ вигляді.                         Корявський В. А.
 *
 *  Ділянка між вузлами A і B задається характеристичним трикутником
 *  A–C–B (C — точка перетину дотичних у кінцях) та інженерним
 *  дискримінантом f. Плечова точка D лежить на медіані CM трикутника:
 *        D = M + f·(C − M),   M = (A + B) / 2,
 *  а сама дуга — раціональна квадратична крива з вагою w = f / (1 − f):
 *        r(t) = [(1−t)²A + 2w·t(1−t)·C + t²B] / [(1−t)² + 2w·t(1−t) + t²].
 *  f < ½ — дуга еліпса, f = ½ — параболи, ½ < f < 1 — гіперболи.
 *  Центр кривої лежить на прямій CM:  O = M − f²/(1 − 2f) · (C − M).
 *
 *  Архітектура (ООП, як у лаб. 1–2): Vec, Mat3, Quadrature, ConicArc,
 *  Contour, Morph, Scene (SVG) та об'єкт-застосунок app.
 * ==================================================================== */

const PI = Math.PI;
const rad = d => d * PI / 180, deg = r => r * 180 / PI;
const $ = id => document.getElementById(id);
const NS = 'http://www.w3.org/2000/svg';
const clamp = (v, lo, hi) => Math.min(hi, Math.max(lo, v));

/* ------------------------------------------------------------- Vec */
class Vec {
    static add(a, b) { return [a[0] + b[0], a[1] + b[1]]; }
    static sub(a, b) { return [a[0] - b[0], a[1] - b[1]]; }
    static mul(a, k) { return [a[0] * k, a[1] * k]; }
    static lerp(a, b, t) { return [a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t]; }
    static mid(a, b) { return [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2]; }
    static dot(a, b) { return a[0] * b[0] + a[1] * b[1]; }
    static cross(a, b) { return a[0] * b[1] - a[1] * b[0]; }
    static len(a) { return Math.hypot(a[0], a[1]); }
    static unit(a) { const l = Vec.len(a) || 1; return [a[0] / l, a[1] / l]; }
}

/* ------------------------------------------------------------ Mat3 */
class Mat3 {
    constructor(m = [1, 0, 0, 0, 1, 0, 0, 0, 1]) { this.m = m; }
    static translate(dx, dy) { return new Mat3([1, 0, dx, 0, 1, dy, 0, 0, 1]); }
    static rotate(d, cx, cy) {
        const c = Math.cos(rad(d)), s = Math.sin(rad(d));
        return new Mat3([c, -s, cx - c * cx + s * cy, s, c, cy - s * cx - c * cy, 0, 0, 1]);
    }
    /** осьова симетрія відносно прямої через (cx, cy) з кутом нахилу ψ */
    static reflect(psi, cx, cy) {
        const c = Math.cos(rad(2 * psi)), s = Math.sin(rad(2 * psi));
        return new Mat3([c, s, cx - c * cx - s * cy, s, -c, cy - s * cx + c * cy, 0, 0, 1]);
    }
    /** this ∘ o (спочатку o) */
    mul(o) {
        const a = this.m, b = o.m, r = new Array(9);
        for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) r[i * 3 + j] = a[i * 3] * b[j] + a[i * 3 + 1] * b[3 + j] + a[i * 3 + 2] * b[6 + j];
        return new Mat3(r);
    }
    apply([x, y]) { const m = this.m; return [m[0] * x + m[1] * y + m[2], m[3] * x + m[4] * y + m[5]]; }
    invert() {
        const [a, b, tx, c, d, ty] = this.m, det = a * d - b * c;
        const ia = d / det, ib = -b / det, ic = -c / det, id = a / det;
        return new Mat3([ia, ib, -(ia * tx + ib * ty), ic, id, -(ic * tx + id * ty), 0, 0, 1]);
    }
    get det() { return this.m[0] * this.m[4] - this.m[1] * this.m[3]; }
}

/* ------------------------------------------------- чисельні методи */
class Quadrature {
    static X = [0, 0.5384693101056831, -0.5384693101056831, 0.9061798459386640, -0.9061798459386640];
    static Wt = [0.5688888888888889, 0.4786286704993665, 0.4786286704993665, 0.2369268850561891, 0.2369268850561891];
    /** складена 5-вузлова квадратура Гауса–Лежандра */
    static gauss(f, lo, hi, panels = 24) {
        const h = (hi - lo) / panels; let sum = 0;
        for (let p = 0; p < panels; p++) {
            const mid = lo + (p + 0.5) * h;
            for (let k = 0; k < 5; k++) sum += Quadrature.Wt[k] * f(mid + Quadrature.X[k] * h / 2);
        }
        return sum * h / 2;
    }
}

/* ------------------------------------------- дуга 2-го порядку */
class ConicArc {
    static NAMES = { el: 'еліпс', pa: 'парабола', hy: 'гіпербола' };
    constructor(A, C, B, f) { this.A = A; this.C = C; this.B = B; this.f = f; }
    get w() { return this.f / (1 - this.f); }
    get M() { return Vec.mid(this.A, this.B); }
    /** плечова точка D на медіані CM */
    get D() { return Vec.lerp(this.M, this.C, this.f); }
    get kind() { return Math.abs(this.f - 0.5) < 0.004 ? 'pa' : this.f < 0.5 ? 'el' : 'hy'; }

    point(t) {
        const s = 1 - t, b0 = s * s, b1 = 2 * this.w * t * s, b2 = t * t, W = b0 + b1 + b2;
        return [(b0 * this.A[0] + b1 * this.C[0] + b2 * this.B[0]) / W, (b0 * this.A[1] + b1 * this.C[1] + b2 * this.B[1]) / W];
    }
    /** r′(t) = (N′·W − N·W′) / W² */
    deriv(t) {
        const s = 1 - t, w = this.w;
        const b0 = s * s, b1 = 2 * w * t * s, b2 = t * t, W = b0 + b1 + b2;
        const d0 = -2 * s, d1 = 2 * w * (1 - 2 * t), d2 = 2 * t, dW = d0 + d1 + d2;
        const out = [0, 0];
        for (const k of [0, 1]) {
            const N = b0 * this.A[k] + b1 * this.C[k] + b2 * this.B[k], dN = d0 * this.A[k] + d1 * this.C[k] + d2 * this.B[k];
            out[k] = (dN * W - N * dW) / (W * W);
        }
        return out;
    }
    sample(n = 24) { return Array.from({ length: n + 1 }, (_, i) => this.point(i / n)); }
    length() { return Quadrature.gauss(t => Vec.len(this.deriv(t)), 0, 1); }
    /** площа сегмента між дугою та хордою AB: ½∮(x dy − y dx) */
    segmentArea() {
        const arc = Quadrature.gauss(t => { const p = this.point(t), d = this.deriv(t); return Vec.cross(p, d); }, 0, 1) / 2;
        return Math.abs(arc + Vec.cross(this.B, this.A) / 2);
    }
    center() {
        if (this.kind === 'pa') return null;
        const k = this.f * this.f / (1 - 2 * this.f);
        return Vec.sub(this.M, Vec.mul(Vec.sub(this.C, this.M), k));
    }
    /** спряжені півдіаметри e1 (уздовж медіани) та e2 (паралельний хорді) → півосі за Аполлонієм */
    axes() {
        const O = this.center(); if (!O) return null;
        const e1 = Vec.sub(this.D, O), c = this.f / (1 - this.f), hB = Vec.sub(this.B, this.M);
        if (this.kind === 'el') {
            const e2 = Vec.mul(hB, 1 / Math.sqrt(1 - c * c));
            const p = Vec.dot(e1, e1) + Vec.dot(e2, e2), q = Math.abs(Vec.cross(e1, e2)), r = Math.sqrt(Math.max(0, p * p - 4 * q * q));
            const a = Math.sqrt((p + r) / 2), b = Math.sqrt(Math.max(0, (p - r) / 2));
            return { O, a, b, e: Math.sqrt(1 - (b * b) / (a * a)) };
        }
        const e2 = Vec.mul(hB, 1 / Math.sqrt(c * c - 1));
        const d = Vec.dot(e1, e1) - Vec.dot(e2, e2), q = Math.abs(Vec.cross(e1, e2));
        const a2 = (d + Math.sqrt(d * d + 4 * q * q)) / 2, a = Math.sqrt(a2), b = Math.sqrt(Math.max(0, a2 - d));
        return { O, a, b, e: Math.sqrt(1 + (b * b) / (a * a)) };
    }
}

/* ---------------------------------------------------- контур */
class Contour {
    /** part: { key, title, closed, knots: [[x, y]], arcs: [{ C: [x, y], f }] } */
    constructor(parts) { this.parts = parts; }
    static fromTable(tbl) {
        return new Contour(Object.entries(tbl).map(([key, p]) => ({
            key, title: p.title, closed: p.closed,
            knots: p.knots.map(q => [...q]), arcs: p.arcs.map(([cx, cy, f]) => ({ C: [cx, cy], f }))
        })));
    }
    static fromJSON(obj) { return new Contour(obj.parts.map(p => ({ ...p, knots: p.knots.map(q => [...q]), arcs: p.arcs.map(a => ({ C: [...a.C], f: a.f })) }))); }
    clone() { return Contour.fromJSON(this); }
    nextKnot(pi, i) { const p = this.parts[pi]; return (i + 1) % p.knots.length; }
    arc(pi, i) { const p = this.parts[pi]; return new ConicArc(p.knots[i], p.arcs[i].C, p.knots[this.nextKnot(pi, i)], p.arcs[i].f); }
    *each() { for (let pi = 0; pi < this.parts.length; pi++) for (let i = 0; i < this.parts[pi].arcs.length; i++) yield { pi, i, part: this.parts[pi], arc: this.arc(pi, i) }; }
    get arcCount() { return this.parts.reduce((s, p) => s + p.arcs.length, 0); }
    /** покрокова лінійна інтерполяція всіх параметрів двох контурів однакової топології */
    static mix(a, b, e) {
        return new Contour(a.parts.map((p, k) => ({
            ...p,
            knots: p.knots.map((q, i) => Vec.lerp(q, b.parts[k].knots[i], e)),
            arcs: p.arcs.map((r, i) => ({ C: Vec.lerp(r.C, b.parts[k].arcs[i].C, e), f: r.f + (b.parts[k].arcs[i].f - r.f) * e }))
        })));
    }
    toJSON() { return { parts: this.parts }; }
}

/* ------------------------------------------- цільові фігури анімації */
class Morph {
    /** щільна замкнена ламана фігури одиничного розміру з центром у (0; 0) */
    static shape(kind) {
        const N = 720, out = [];
        for (let i = 0; i < N; i++) {
            const t = 2 * PI * i / N;
            if (kind === 'hearts') out.push([16 * Math.sin(t) ** 3 / 17, (13 * Math.cos(t) - 5 * Math.cos(2 * t) - 2 * Math.cos(3 * t) - Math.cos(4 * t) + 3) / 17]);
            else if (kind === 'circles') out.push([Math.cos(t), Math.sin(t)]);
            else {                                                        // рівносторонній трикутник
                const k = Math.floor(3 * i / N), u = 3 * i / N - k, a0 = rad(90 + 120 * k), a1 = rad(90 + 120 * (k + 1));
                out.push([Math.cos(a0) + (Math.cos(a1) - Math.cos(a0)) * u, Math.sin(a0) + (Math.sin(a1) - Math.sin(a0)) * u]);
            }
        }
        return out;
    }
    /** точка ламаної на промені з центру під полярним кутом φ */
    static onRay(poly, phi) {
        const d = [Math.cos(phi), Math.sin(phi)];
        let best = null;
        for (let i = 0; i < poly.length; i++) {
            const p = poly[i], q = poly[(i + 1) % poly.length], e = Vec.sub(q, p), den = Vec.cross(d, e);
            if (Math.abs(den) < 1e-12) continue;
            const t = Vec.cross(p, e) / den, u = Vec.cross(p, d) / den;
            if (t > 0 && u >= 0 && u <= 1 && (!best || t > best.t)) best = { t, i, u };
        }
        return best;
    }
    /** рівномірний за довжиною розподіл n точок по відкритій ламаній */
    static along(poly, n) {
        const cum = [0]; for (let i = 1; i < poly.length; i++) cum.push(cum[i - 1] + Vec.len(Vec.sub(poly[i], poly[i - 1])));
        const L = cum[cum.length - 1], out = []; let k = 1;
        for (let j = 0; j < n; j++) {
            const s = L * j / (n - 1); while (k < poly.length - 1 && cum[k] < s) k++;
            const t = (s - cum[k - 1]) / (cum[k] - cum[k - 1] || 1); out.push({ p: Vec.lerp(poly[k - 1], poly[k], t), i: k - 1 });
        }
        return out;
    }
    /** вершина C і дискримінант f дуги між точками ламаної з відомими дотичними */
    static arcFor(A, tA, B, tB) {
        const L = Vec.len(Vec.sub(B, A)), den = Vec.cross(tA, tB);
        const turn = Math.acos(clamp(Vec.dot(tA, tB), -1, 1)), h = turn / 2;
        if (Math.abs(den) > 1e-9 && turn > rad(1.5)) {
            const r = Vec.sub(B, A), s = Vec.cross(r, tB) / den, q = -Vec.cross(r, tA) / den;   // C = A + s·tA = B − q·tB
            if (s > 0 && q > 0 && s < 3 * L && q < 3 * L) return { C: Vec.add(A, Vec.mul(tA, s)), f: Math.cos(h) / (1 + Math.cos(h)) };
        }
        return { C: Vec.mid(A, B), f: 0.5 };                                // відрізок прямої
    }
    static build(src, kind) {
        const base = Morph.shape(kind);
        return new Contour(src.parts.map(part => {
            const xs = part.knots.map(q => q[0]), ys = part.knots.map(q => q[1]);
            const cx = (Math.min(...xs) + Math.max(...xs)) / 2, cy = (Math.min(...ys) + Math.max(...ys)) / 2;
            const R = Math.max(Math.max(...xs) - Math.min(...xs), Math.max(...ys) - Math.min(...ys)) * (part.closed ? 0.5 : 0.62);
            const poly = base.map(([x, y]) => [cx + R * x, cy + R * y]);
            const tang = i => Vec.unit(Vec.sub(poly[(i + 1) % poly.length], poly[(i - 1 + poly.length) % poly.length]));
            const n = part.knots.length;
            let pts;
            if (part.closed) {
                // старт — полярний кут першого вузла, напрям обходу — як у вихідної ділянки
                const area = part.knots.reduce((s, q, i) => s + Vec.cross(q, part.knots[(i + 1) % n]), 0), dir = Math.sign(area) || 1;
                const phi0 = Math.atan2(part.knots[0][1] - cy, part.knots[0][0] - cx);
                const local = base.map(q => Vec.mul(q, R));
                pts = Array.from({ length: n }, (_, j) => {
                    const hit = Morph.onRay(local, phi0 + dir * 2 * PI * j / n);
                    return { p: [cx + hit.t * Math.cos(phi0 + dir * 2 * PI * j / n), cy + hit.t * Math.sin(phi0 + dir * 2 * PI * j / n)], i: hit.i };
                });
                pts.forEach(o => o.t = tang(o.i));
            } else {
                // відкрита ділянка (шнур) — фігура з розривом угорі
                const k0 = Math.round(base.length * 0.04), open = poly.slice(k0, base.length - k0);
                pts = Morph.along(open, n);
                pts.forEach(o => o.t = tang(o.i + k0));
            }
            // дотичні орієнтуємо за напрямом обходу нових вузлів
            pts.forEach((o, j) => {
                const nx = pts[Math.min(j + 1, n - 1)].p, pv = pts[Math.max(j - 1, 0)].p;
                const chord = part.closed ? Vec.sub(pts[(j + 1) % n].p, pts[(j - 1 + n) % n].p) : Vec.sub(nx, pv);
                if (Vec.dot(o.t, chord) < 0) o.t = Vec.mul(o.t, -1);
            });
            const arcs = part.arcs.map((_, i) => { const a = pts[i], b = pts[(i + 1) % n]; return Morph.arcFor(a.p, a.t, b.p, b.t); });
            return { ...part, knots: pts.map(o => o.p), arcs };
        }));
    }
}

/* ------------------------------------------------------------ стан */
const state = {
    contour: Contour.fromTable(ROPE_TABLE),
    sc: 30, view: { ox: 0, oy: 0 },
    show: { poly: true, knots: true, shoulder: true, kind: false, grid: true },
    sel: { pi: 1, i: 6 },
    shift: [0, 0], rot: { x: 0, y: 0, a: 0 }, refl: { on: false, psi: 90 },
    morph: { tau: 0, playing: false, dir: 1, last: 0, target: null, kind: 'hearts' },
    dragging: false
};
const eucl = () => {
    let T = Mat3.rotate(state.rot.a, state.rot.x, state.rot.y).mul(Mat3.translate(state.shift[0], state.shift[1]));
    if (state.refl.on) T = Mat3.reflect(state.refl.psi, state.rot.x, state.rot.y).mul(T);
    return T;
};

/* ---------------------------------------------------- SVG-сцена */
class Scene {
    constructor() {
        this.svg = $('scene');
        this.gridG = this.mk('g', {}, this.svg);
        this.world = this.mk('g', {}, this.svg);
        this.ghostG = this.mk('g', {}, this.world);
        this.polyG = this.mk('g', {}, this.world);
        this.curveG = this.mk('g', {}, this.world);
        this.ui = this.mk('g', {}, this.svg);
        this.firstFit = true;
        new ResizeObserver(() => this.fit()).observe($('svgWrap'));
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
        if (this.firstFit && this.W > 0) { this.firstFit = false; app.fitView(); }
        app.draw();
    }
    toScreen([x, y]) { return [state.view.ox + x * state.sc, state.view.oy - y * state.sc]; }
    toWorld(sx, sy) { return [(sx - state.view.ox) / state.sc, (state.view.oy - sy) / state.sc]; }
    d(pts) { return pts.map((p, i) => (i ? 'L' : 'M') + p[0].toFixed(4) + ' ' + p[1].toFixed(4)).join(''); }
    text(s, x, y, cls) { this.mk('text', { x, y, class: cls }, this.ui).textContent = s; }

    drawGrid() {
        this.gridG.replaceChildren();
        if (!state.show.grid) return;
        const [x0, y1] = this.toWorld(0, 0), [x1, y0] = this.toWorld(this.W, this.H), line = (a, b, cls) => {
            const p = this.toScreen(a), q = this.toScreen(b); this.mk('line', { x1: p[0], y1: p[1], x2: q[0], y2: q[1], class: cls }, this.gridG);
        };
        for (let i = Math.floor(x0); i <= Math.ceil(x1); i++) line([i, y0], [i, y1], i % 5 ? 'gridl' : 'gridM');
        for (let j = Math.floor(y0); j <= Math.ceil(y1); j++) line([x0, j], [x1, j], j % 5 ? 'gridl' : 'gridM');
        line([x0, 0], [x1, 0], 'axisX'); line([0, y0], [0, y1], 'axisY');
        const o = this.toScreen([0, 0]), lbl = (s, x, y) => this.mk('text', { x, y, class: 'lbl' }, this.gridG).textContent = s;
        for (let i = Math.ceil(x0); i <= x1; i++) if (i && i % 2 === 0) lbl(i, o[0] + i * state.sc - 4, o[1] + 15);
        for (let j = Math.ceil(y0); j <= y1; j++) if (j && j % 2 === 0) lbl(j, o[0] + 7, o[1] - j * state.sc + 4);
    }
}

/* ------------------------------------------------------------- App */
const f3 = v => (Math.round(v * 1e3) / 1e3).toFixed(3);
const f4 = v => (Math.round(v * 1e4) / 1e4).toString();
const pt = p => `(${f3(p[0])}; ${f3(p[1])})`;

const app = {
    /** контур, що зараз показано: вихідний або проміжний кадр переходу */
    shown() {
        const tau = state.morph.tau;
        if (tau <= 0) return state.contour;
        const e = 0.5 - 0.5 * Math.cos(PI * tau);               // косинусне згладжування ходу
        return Contour.mix(state.contour, this.target(), e);
    },
    target() { return state.morph.target || (state.morph.target = Morph.build(state.contour, state.morph.kind)); },
    editable() { return state.morph.tau === 0; },

    draw() {
        const sc = scene; if (!sc.W) return;
        const T = eucl(), shown = this.shown();
        sc.world.setAttribute('transform', `matrix(${state.sc} 0 0 ${-state.sc} ${state.view.ox} ${state.view.oy})`);
        sc.drawGrid();
        for (const g of [sc.ghostG, sc.polyG, sc.curveG, sc.ui]) g.replaceChildren();

        // цільова фігура
        if ($('ghost').checked) for (const { arc } of this.target().each()) sc.mk('path', { d: sc.d(arc.sample(16).map(p => T.apply(p))), class: 'ghostArc' }, sc.ghostG);

        // дуги контуру + невидимі «товсті» копії для вибору кліком
        for (const { pi, i, arc } of shown.each()) {
            const d = sc.d(arc.sample(24).map(p => T.apply(p))), sel = pi === state.sel.pi && i === state.sel.i;
            sc.mk('path', { d, class: 'arc' + (state.show.kind ? ' ' + arc.kind : '') + (sel ? ' sel' : '') }, sc.curveG);
            sc.mk('path', { d, class: 'hitArc', 'data-kind': 'arc', 'data-pi': pi, 'data-i': i }, sc.curveG);
            if (state.show.poly) sc.mk('polyline', { points: [arc.A, arc.C, arc.B].map(p => T.apply(p).join(',')).join(' '), class: 'poly' + (sel ? ' sel' : '') }, sc.polyG);
        }

        // вузли та вершини C
        if (state.show.knots) shown.parts.forEach((part, pi) => {
            if (state.show.poly) part.arcs.forEach((a, i) => { const q = sc.toScreen(T.apply(a.C)); sc.mk('circle', { cx: q[0], cy: q[1], r: pi === state.sel.pi && i === state.sel.i ? 5 : 3.6, class: 'pole', 'data-kind': 'pole', 'data-pi': pi, 'data-i': i }, sc.ui); });
            part.knots.forEach((k, i) => { const q = sc.toScreen(T.apply(k)); sc.mk('circle', { cx: q[0], cy: q[1], r: 4, class: 'knot', 'data-kind': 'knot', 'data-pi': pi, 'data-i': i }, sc.ui); });
        });

        // вибрана дуга: медіана CM і плечова точка D (поверх усіх маркерів, щоб її можна було схопити)
        const cur = shown.parts[state.sel.pi] ? shown.arc(state.sel.pi, state.sel.i) : null;
        if (cur && state.show.shoulder && state.show.poly) {
            const C = sc.toScreen(T.apply(cur.C)), M = sc.toScreen(T.apply(cur.M)), D = sc.toScreen(T.apply(cur.D));
            sc.mk('line', { x1: C[0], y1: C[1], x2: M[0], y2: M[1], class: 'median' }, sc.ui);
            const A = sc.toScreen(T.apply(cur.A)), B = sc.toScreen(T.apply(cur.B));
            sc.mk('line', { x1: A[0], y1: A[1], x2: B[0], y2: B[1], class: 'median', opacity: 0.45 }, sc.ui);
            sc.mk('rect', { x: D[0] - 5, y: D[1] - 5, width: 10, height: 10, transform: `rotate(45 ${D[0]} ${D[1]})`, class: 'shoulder', 'data-kind': 'shoulder' }, sc.ui);
            for (const [s, p] of [['A', A], ['B', B], ['C', C], ['D', D], ['M', M]]) sc.text(s, p[0] + 8, p[1] - 8, 'tag');
        }

        // центр обертання та вісь симетрії
        if (state.rot.a || state.refl.on) {
            const q = sc.toScreen([state.rot.x, state.rot.y]);
            sc.mk('path', { d: `M${q[0] - 10} ${q[1]}H${q[0] + 10}M${q[0]} ${q[1] - 10}V${q[1] + 10}`, class: 'pivot' }, sc.ui);
            sc.mk('circle', { cx: q[0], cy: q[1], r: 6, class: 'pivot' }, sc.ui);
            sc.text('P', q[0] + 9, q[1] + 16, 'mk');
            if (state.refl.on) {
                const d = [Math.cos(rad(state.refl.psi)), Math.sin(rad(state.refl.psi))], L = 1e3;
                const a = sc.toScreen([state.rot.x - d[0] * L, state.rot.y - d[1] * L]), b = sc.toScreen([state.rot.x + d[0] * L, state.rot.y + d[1] * L]);
                sc.mk('line', { x1: a[0], y1: a[1], x2: b[0], y2: b[1], class: 'reflAxis' }, sc.ui);
            }
        }
        this.panels(shown);
    },

    /* ----- нижня панель та бокові підписи ----- */
    panels(shown) {
        const part = state.contour.parts[state.sel.pi], arc = shown.arc(state.sel.pi, state.sel.i);
        $('selCap').textContent = `${part.title}: дуга ${state.sel.i + 1} з ${part.arcs.length} · ${ConicArc.NAMES[arc.kind]}`;
        if (document.activeElement !== $('f')) $('f').value = +arc.f.toFixed(3);
        $('fR').value = arc.f;
        if (!$('dTbl').hidden && !state.dragging) this.fillTable(shown);
        if (!$('dArc').hidden) this.fillArcCards(shown, arc);
        $('pOut').textContent = state.morph.tau.toFixed(3);
    },
    fillTable(shown) {
        let n = 0; const rows = [];
        for (const { pi, i, part, arc } of shown.each()) {
            n++;
            const sel = pi === state.sel.pi && i === state.sel.i ? ' class="sel"' : '';
            rows.push(`<tr${sel} data-pi="${pi}" data-i="${i}"><td>${n}</td><td>${part.title} · ${i + 1}</td><td>${pt(arc.A)}</td><td>${pt(arc.C)}</td><td>${pt(arc.B)}</td><td>${arc.f.toFixed(3)}</td><td>${arc.w.toFixed(3)}</td><td class="k-${arc.kind}">${ConicArc.NAMES[arc.kind]}</td><td>${f3(arc.length())}</td></tr>`);
        }
        $('tbl').tBodies[0].innerHTML = rows.join('');
        const tr = $('tbl').tBodies[0].querySelector('tr.sel');
        if (tr && this.scrollToSel) { tr.scrollIntoView({ block: 'nearest' }); this.scrollToSel = false; }
    },
    fillArcCards(shown, arc) {
        const part = shown.parts[state.sel.pi], ax = arc.axes(), T = eucl(), O = arc.center();
        const card = (title, big, small = '') => `<div class="stat"><small>${title}</small><b>${big}</b>${small ? `<p>${small}</p>` : ''}</div>`;
        const pr = (k, s) => `${f4(arc.A[k])}·(1−t)² ${s} 2·${f4(arc.w)}·${f4(arc.C[k])}·t(1−t) + ${f4(arc.B[k])}·t²`;
        $('dArc').innerHTML = [
            card('Дуга', `${part.title} · ${state.sel.i + 1} / ${part.arcs.length}`, `A = ${pt(arc.A)} → B = ${pt(arc.B)}`),
            card('Тип кривої', ConicArc.NAMES[arc.kind], `f = ${f4(arc.f)}, w = f/(1−f) = ${f4(arc.w)}`),
            card('Вершина C та плечова точка D', pt(arc.C), `D = M + f(C − M) = ${pt(arc.D)}; M = ${pt(arc.M)}`),
            card('Центр кривої O', O ? pt(O) : '∞', O ? `O = M − f²/(1−2f)·(C − M)` : 'у параболи центр — нескінченно віддалена точка осі ∥ CM'),
            card(arc.kind === 'pa' ? 'Ексцентриситет' : 'Півосі та ексцентриситет', ax ? `a = ${f4(ax.a)}, b = ${f4(ax.b)}` : 'e = 1', ax ? `e = ${f4(ax.e)} (${arc.kind === 'el' ? 'e < 1' : 'e > 1'})` : 'вісь параболи паралельна медіані CM'),
            card('Довжина дуги', f4(arc.length()), 's = ∫₀¹ |r′(t)| dt, Гаусс–Лежандр 5×24'),
            card('Площа сегмента «дуга – хорда»', f4(arc.segmentArea()), '½∮(x dy − y dx) по дузі та хорді BA'),
            card('Після перетворення T', `A′ = ${pt(T.apply(arc.A))}`, `B′ = ${pt(T.apply(arc.B))}, C′ = ${pt(T.apply(arc.C))}`),
            `<div class="stat" style="grid-column: span 2"><small>Параметричні рівняння дуги, t ∈ [0; 1]</small><p style="color:var(--ink);font:12px Consolas,monospace">x(t) = [${pr(0, '+')}] / Φ(t)<br>y(t) = [${pr(1, '+')}] / Φ(t)<br>Φ(t) = (1−t)² + 2·${f4(arc.w)}·t(1−t) + t²</p></div>`
        ].join('');
    },

    select(pi, i, scroll = true) { state.sel = { pi, i }; this.scrollToSel = scroll; this.draw(); },
    step(d) {
        const list = [...state.contour.each()].map(o => [o.pi, o.i]);
        let k = list.findIndex(([a, b]) => a === state.sel.pi && b === state.sel.i);
        k = (k + d + list.length) % list.length; this.select(...list[k]);
    },
    setF(f, all = false) {
        f = clamp(f, 0.02, 0.98);
        const part = state.contour.parts[state.sel.pi];
        (all ? part.arcs : [part.arcs[state.sel.i]]).forEach(a => a.f = f);
        state.morph.target = null; this.draw();
    },
    fitView() {
        const T = eucl(), pts = this.shown().parts.flatMap(p => [...p.knots, ...p.arcs.map(a => a.C)].map(q => T.apply(q)));
        const xs = pts.map(q => q[0]), ys = pts.map(q => q[1]);
        const x0 = Math.min(...xs), x1 = Math.max(...xs), y0 = Math.min(...ys), y1 = Math.max(...ys);
        state.sc = clamp(Math.min((scene.W - 80) / (x1 - x0), (scene.H - 60) / (y1 - y0)), 8, 120);
        state.view.ox = scene.W / 2 - (x0 + x1) / 2 * state.sc; state.view.oy = scene.H / 2 + (y0 + y1) / 2 * state.sc;
    },

    /* ----- анімація переходу ----- */
    tick(now) {
        const m = state.morph; if (!m.playing) return;
        const dt = (now - m.last) / 1000; m.last = now;
        m.tau += m.dir * dt / Math.max(0.5, parseFloat($('dur').value) || 4);
        if (m.tau >= 1) { m.tau = 1; if ($('pingpong').checked) m.dir = -1; else this.play(false); }
        if (m.tau <= 0) { m.tau = 0; m.dir = 1; }
        $('pR').value = m.tau; this.draw();
        requestAnimationFrame(n => this.tick(n));
    },
    play(on) {
        const m = state.morph; m.playing = on;
        $('btnPlay').textContent = on ? '❚❚ Пауза' : '▶ Відтворити';
        if (on) { if (m.tau >= 1 && !$('pingpong').checked) m.tau = 0; m.last = performance.now(); requestAnimationFrame(n => this.tick(n)); }
    },

    /* ----- миша: перетягування вузлів, вершин C та плечових точок ----- */
    pointer() {
        const svg = scene.svg; let drag = null;
        const local = e => { const r = svg.getBoundingClientRect(); return [e.clientX - r.left, e.clientY - r.top]; };
        const model = e => eucl().invert().apply(scene.toWorld(...local(e)));
        svg.addEventListener('pointerdown', e => {
            const [sx, sy] = local(e), kind = e.target.dataset.kind;
            if (e.shiftKey) { const w = scene.toWorld(sx, sy); $('rx').value = +w[0].toFixed(2); $('ry').value = +w[1].toFixed(2); this.readTransform(); return; }
            if (kind && kind !== 'arc' && this.editable()) {
                drag = { kind, pi: +e.target.dataset.pi, i: +e.target.dataset.i };
                if (kind !== 'shoulder' && kind !== 'knot') state.sel = { pi: drag.pi, i: drag.i };
                if (kind === 'knot') { const p = state.contour.parts[drag.pi]; state.sel = { pi: drag.pi, i: Math.min(drag.i, p.arcs.length - 1) }; }
                state.dragging = true;
            } else drag = { pan: true, sx, sy, ox: state.view.ox, oy: state.view.oy, moved: 0, target: e.target };
            svg.setPointerCapture(e.pointerId); svg.classList.add('drag');
        });
        svg.addEventListener('pointermove', e => {
            const [sx, sy] = local(e), w = scene.toWorld(sx, sy);
            $('readout').textContent = `x = ${f3(w[0])}   y = ${f3(w[1])}`;
            if (!drag) return;
            if (drag.pan) {
                drag.moved = Math.max(drag.moved, Math.hypot(sx - drag.sx, sy - drag.sy));
                state.view.ox = drag.ox + sx - drag.sx; state.view.oy = drag.oy + sy - drag.sy; this.draw(); return;
            }
            const p = model(e), part = state.contour.parts[drag.pi];
            if (drag.kind === 'knot') part.knots[drag.i] = [+p[0].toFixed(3), +p[1].toFixed(3)];
            else if (drag.kind === 'pole') part.arcs[drag.i].C = [+p[0].toFixed(3), +p[1].toFixed(3)];
            else {                                    // проєкція курсора на медіану CM → новий дискримінант f
                const arc = state.contour.arc(state.sel.pi, state.sel.i), mc = Vec.sub(arc.C, arc.M);
                state.contour.parts[state.sel.pi].arcs[state.sel.i].f = +clamp(Vec.dot(Vec.sub(p, arc.M), mc) / Vec.dot(mc, mc), 0.02, 0.98).toFixed(3);
            }
            state.morph.target = null;
            $('readout').textContent = drag.kind === 'shoulder' ? `f = ${state.contour.parts[state.sel.pi].arcs[state.sel.i].f}` : `${drag.kind === 'knot' ? 'вузол A' : 'вершина C'} → ${pt(p)}`;
            this.draw();
        });
        svg.addEventListener('pointerup', () => {
            if (drag && drag.pan && drag.moved < 4 && drag.target.dataset.kind === 'arc') this.select(+drag.target.dataset.pi, +drag.target.dataset.i);
            drag = null; state.dragging = false; svg.classList.remove('drag'); this.draw();
        });
        svg.addEventListener('wheel', e => {
            e.preventDefault();
            const [sx, sy] = local(e), k = e.deltaY < 0 ? 1.1 : 1 / 1.1, sc = clamp(state.sc * k, 8, 160), r = sc / state.sc;
            state.view.ox = sx - (sx - state.view.ox) * r; state.view.oy = sy - (sy - state.view.oy) * r; state.sc = sc; this.draw();
        }, { passive: false });
        svg.addEventListener('pointerleave', () => { $('readout').textContent = ''; });
    },

    readTransform() {
        const num = id => { const v = parseFloat($(id).value); return Number.isFinite(v) ? v : 0; };
        state.shift = [num('shx'), num('shy')]; state.rot = { x: num('rx'), y: num('ry'), a: num('ra') };
        state.refl = { on: $('refl').checked, psi: num('rpsi') };
        $('raR').value = state.rot.a; $('rpsiR').value = state.refl.psi; this.draw();
    },
    bind() {
        for (const [id, key] of [['vPoly', 'poly'], ['vKnots', 'knots'], ['vShoulder', 'shoulder'], ['vKind', 'kind'], ['vGrid', 'grid']])
            $(id).onchange = e => { state.show[key] = e.target.checked; this.draw(); };
        $('f').oninput = () => { const v = parseFloat($('f').value); if (Number.isFinite(v)) this.setF(v); };
        $('fR').oninput = () => this.setF(+$('fR').value);
        document.querySelectorAll('.trio button').forEach(b => b.onclick = () => this.setF(+b.dataset.f));
        $('btnAllF').onclick = () => this.setF(state.contour.parts[state.sel.pi].arcs[state.sel.i].f, true);
        $('sPrev').onclick = () => this.step(-1); $('sNext').onclick = () => this.step(1);

        for (const id of ['shx', 'shy', 'rx', 'ry', 'ra', 'rpsi']) $(id).oninput = () => this.readTransform();
        $('refl').onchange = () => this.readTransform();
        $('raR').oninput = () => { $('ra').value = $('raR').value; this.readTransform(); };
        $('rpsiR').oninput = () => { $('rpsi').value = $('rpsiR').value; this.readTransform(); };
        $('btnReset').onclick = () => { for (const id of ['shx', 'shy', 'rx', 'ry', 'ra']) $(id).value = 0; $('rpsi').value = 90; $('refl').checked = false; this.readTransform(); };

        $('btnDefault').onclick = () => { state.contour = Contour.fromTable(ROPE_TABLE); state.morph.target = null; this.draw(); };
        $('btnFit').onclick = () => { this.fitView(); this.draw(); };
        $('btnSave').onclick = () => {
            const a = document.createElement('a');
            a.href = URL.createObjectURL(new Blob([JSON.stringify({ method: 'conic-engineering-form', variant: 4, ...state.contour.toJSON() }, null, 1)], { type: 'application/json' }));
            a.download = 'rope_contour.json'; a.click(); URL.revokeObjectURL(a.href);
        };
        $('fileIn').onchange = async e => {
            const file = e.target.files[0]; if (!file) return;
            try { state.contour = Contour.fromJSON(JSON.parse(await file.text())); state.sel = { pi: 0, i: 0 }; state.morph.target = null; this.fitView(); this.draw(); }
            catch (err) { alert('Файл не є контуром цієї програми: ' + err.message); }
            e.target.value = '';
        };

        document.querySelectorAll('#dockTabs button').forEach(b => b.onclick = () => {
            document.querySelectorAll('#dockTabs button').forEach(x => x.classList.toggle('on', x === b));
            $('dTbl').hidden = b.dataset.d !== 'tbl'; $('dArc').hidden = b.dataset.d !== 'arc'; $('dAnim').hidden = b.dataset.d !== 'anim';
            this.scrollToSel = true; this.draw();
        });
        $('tbl').tBodies[0].onclick = e => { const tr = e.target.closest('tr'); if (tr) this.select(+tr.dataset.pi, +tr.dataset.i, false); };

        $('target').onchange = () => { state.morph.kind = $('target').value; state.morph.target = null; this.draw(); };
        $('pR').oninput = () => { state.morph.tau = +$('pR').value; this.draw(); };
        $('btnPlay').onclick = () => this.play(!state.morph.playing);
        $('btnBack').onclick = () => { this.play(false); state.morph.tau = 0; state.morph.dir = 1; $('pR').value = 0; this.draw(); };
        $('ghost').onchange = () => this.draw();

        window.addEventListener('keydown', e => {
            if (/INPUT|SELECT/.test(e.target.tagName)) return;
            if (e.key === 'p' || e.key === 'P' || e.key === 'з' || e.key === 'З') { $('vPoly').checked = state.show.poly = !state.show.poly; this.draw(); }
            if (e.key === 'ArrowRight') this.step(1);
            if (e.key === 'ArrowLeft') this.step(-1);
        });
        this.pointer();
    },
    applyQuery() {
        const q = new URLSearchParams(location.search).get('s'); if (!q) return;
        try {
            const s = JSON.parse(q), c = state.contour;
            if (s.sel) state.sel = { pi: s.sel[0], i: s.sel[1] };
            if (s.knots) for (const [pi, i, x, y] of s.knots) c.parts[pi].knots[i] = [x, y];
            if (s.poles) for (const [pi, i, x, y] of s.poles) c.parts[pi].arcs[i].C = [x, y];
            if (s.f) for (const [pi, i, f] of s.f) c.parts[pi].arcs[i].f = f;
            if (s.show) for (const k in s.show) { state.show[k] = s.show[k]; const id = { poly: 'vPoly', knots: 'vKnots', shoulder: 'vShoulder', kind: 'vKind', grid: 'vGrid' }[k]; $(id).checked = s.show[k]; }
            if (s.target) { $('target').value = s.target; state.morph.kind = s.target; }
            if (s.tau !== undefined) { state.morph.tau = s.tau; $('pR').value = s.tau; }
            if (s.ghost) $('ghost').checked = true;
            if (s.shift) { $('shx').value = s.shift[0]; $('shy').value = s.shift[1]; }
            if (s.rot) { $('rx').value = s.rot.x; $('ry').value = s.rot.y; $('ra').value = s.rot.a; }
            if (s.refl) { $('refl').checked = true; $('rpsi').value = s.refl; }
            if (s.open) document.querySelectorAll('#rail details').forEach(d => d.open = s.open.includes(d.querySelector('summary').textContent.slice(0, 4)));
            this.readTransform();
            if (s.tab) document.querySelector(`#dockTabs button[data-d="${s.tab}"]`).click();
            if (s.play) this.play(true);
            this.refit = s.refit;
        } catch (e) { /* некоректний запит ігноруємо */ }
    }
};

const scene = new Scene();
app.bind(); app.applyQuery();
scene.fit();
if (app.refit) { app.fitView(); app.draw(); }
