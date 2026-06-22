# ballistX

FLOSS ballistics calculator for competitive shooting events.

## Features

- 3-DOF point-mass trajectory solver.
- Modular ammunition catalog.
- Automatic environmental corrections.
- Wind, material, and construction adjustments.
- Barrel-specific adjustments.
- Imperial/Metric toggle & MOA/MIL toggle.
- Clamping for reasonable ranges.

## Building & running

Requires Go 1.26+ and on Linux-based systems, `libGL`, `libX11`, `libXcursor`, `libXrandr`, `libXinerama`, `libXi`.

```sh
go build -o ballistx .
./ballistx
```

or just:

```sh
go run .
```

## How to Use ballistX

1. Pick ammunition from the dropdown, or click "Manage Ammo Presets" to add your own load (brand, line, cartridge, bullet weight, BC, drag function, muzzle velocity, test barrel length, bullet diameter/length, construction, and a reasonable max range).
2. Set your conditions.
3. Click "Calculate Trajectory". The chart shows drop and wind drift vs. range; the table gives a full row-by-row breakdown with MOA/MIL corrections.
4. Switch units anytime with the Imperial/Metric radio. Already-entered values convert automatically.

Your ammo catalog is stored at `~/.config/ballistx/catalog.json`.

### Modeling notes & limitations

This is a standard point-mass model — the same class of math used by most consumer ballistic calculators (no Coriolis effect, no aerodynamic jump). A few specific approximations worth knowing about:

- Drag tables are the classic public-domain G1/G7 standard projectile tables (Ingalls/McCoy/Litz lineage), not a measurement of any specific bullet. Your BC scales how closely your bullet tracks the chosen standard.
- Barrel-length velocity adjustment is a linear fps-per-inch estimate relative to the preset's test-barrel length. Real velocity-vs-length curves are charge- and cartridge-specific; treat this as a reasonable correction, not gospel.
- Spin drift uses the Litz empirical approximation from gyroscopic stability (Miller twist rule), not a full 6-DOF simulation.
- "Material Adjustment" (bullet construction) and the manual BC% fine-tune are provided as user-facing knobs since jacket/construction effects on BC retention are bullet-specific and not something a generic model can derive — correct them against your own chronograph/drop data where it matters.

The included physics has a regression test (`internal/ballistics/sanity_test.go`) checking solved drop at 500/1000yd for a well-documented .308 Win 168gr load against published reference bands — run `go test ./...` to verify after any changes to the solver.
