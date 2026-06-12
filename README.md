# Correlation Tightness

A small interactive tool for exploring **correlation matrix tightness**: given a
correlation matrix and one off-diagonal entry, how constrained is that entry by
the rest of the matrix?

🔗 **Live site:** <https://correlation-tightness.com>

- **Explore** — generate or upload a correlation matrix, click any off-diagonal
  cell to see its feasible range given the other entries, then drag or type a
  value and watch the other entries' ranges shift in response.
- **Theory** — the derivation behind the bounds: feasible range via the Schur
  complement, and its connection to partial correlation.

All computation runs client-side in vanilla JavaScript — no backend, no network
round-trips. The math primitives (seeded RNG, Wishart correlation sampling,
Cholesky solve, Jacobi eigendecomposition, Schur-complement bounds) live in a
single dependency-free file, `static/js/linalg.js`.

## Built with

- [Zola](https://www.getzola.org/) — static site generator
- [lightspeed](https://github.com/carpetscheme/lightspeed) theme (MIT, vendored under `themes/`)
- [KaTeX](https://katex.org/) 0.16.11 (MIT, vendored under `static/katex/`) — math rendering on the theory page

## Running locally

Requires [Zola](https://www.getzola.org/documentation/getting-started/installation/) (0.22.x).

```sh
zola serve --port 1111   # hot-reloads on file changes
```

Then open <http://127.0.0.1:1111/>.

## Building

```sh
zola build               # output → public/
```

`public/` is a self-contained static bundle, deployable to any static host.

## License

[MIT](./LICENSE) for the project's own code. Vendored dependencies retain their
own licenses: the lightspeed theme (`themes/lightspeed/LICENSE.md`) and KaTeX,
both MIT.
