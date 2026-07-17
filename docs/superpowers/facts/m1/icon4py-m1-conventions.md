# icon4py facts for the M1 plan (amps package foundations)

Worktree: `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics` (branch `amps_microphysics`)

---

## 1. muphys constants pattern (float-Enum)

File: `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/model/atmosphere/subgrid_scale_physics/muphys/src/icon4py/model/atmosphere/subgrid_scale_physics/muphys/core/common/constants.py` — quoted in full:

```python
# ICON4Py - ICON inspired code in Python and GT4Py
#
# Copyright (c) 2022-2024, ETH Zurich and MeteoSwiss
# All rights reserved.
#
# Please, refer to the LICENSE file in the root directory.
# SPDX-License-Identifier: BSD-3-Clause

import enum

from icon4py.model.common import type_alias as ta


class GraupelConsts(ta.wpfloat, enum.Enum):
    rho_00 = 1.225  # reference air density
    q1 = 8.0e-6
    qmin = 1.0e-15  # threshold for computation
    ams = 0.069  # Formfactor in the mass-size relation of snow particles
    bms = 2.0  # Exponent in the mass-size relation of snow particles
    v0s = 25.0  # prefactor in snow fall speed
    v1s = 0.5  # Exponent in the terminal velocity for snow
    m0_ice = 1.0e-12  # initial crystal mass for cloud ice nucleation
    ci = 2108.0  # specific heat of ice
    tx = 3339.5
    tfrz_het1 = 267.15  # temperature for het. freezing of cloud water with supersat => TMELT - 6.0
    tfrz_het2 = 248.15  # temperature for het. freezing of cloud water => TMELT - 25.0
    tfrz_hom = 236.15  # temperature for hom. freezing of cloud water => TMELT - 37.0
    lvc = 3135383.2031928  # invariant part of vaporization enthalpy => alv - (cpv - clw) * tmelt
    lsc = 2899657.201  # invariant part of vaporization enthalpy => als - (cpv - ci) * tmelt


class ThermodynamicConsts(ta.wpfloat, enum.Enum):
    # Thermodynamic constants for the dry and moist atmosphere
    # Dry air
    rd = 287.04  # [J/K/kg] gas constant
    cpd = 1004.64  # [J/K/kg] specific heat at constant pressure
    cvd = 717.60  # [J/K/kg] specific heat at constant volume => cpd - rd
    con_m = 1.50e-5  # [m^2/s]  kinematic viscosity of dry air
    con_h = 2.20e-5  # [m^2/s]  scalar conductivity of dry air
    con0_h = 2.40e-2  # [J/m/s/K] thermal conductivity of dry air
    eta0d = 1.717e-5  # [N*s/m2] dyn viscosity of dry air at tmelt
    # H2O
    # gas
    rv = 461.51  # [J/K/kg] gas constant for water vapor
    cpv = 1869.46  # [J/K/kg] specific heat at constant pressure
    cvv = 1407.95  # [J/K/kg] specific heat at constant volume => cpv - rv
    dv0 = 2.22e-5  # [m^2/s]  diff coeff of H2O vapor in dry air at tmelt
    # liquid / water
    rhoh2o = 1000.0  # [kg/m3]  density of liquid water
    # solid / ice
    rhoice = 916.7  # [kg/m3]  density of pure ice
    cv_i = 2000.0
    # phase changes
    alv = 2.5008e6  # [J/kg]   latent heat for vaporisation
    als = 2.8345e6  # [J/kg]   latent heat for sublimation
    alf = 333700.0  # [J/kg]   latent heat for fusion => als - alv
    tmelt = 273.15  # [K]      melting temperature of ice/snow
    t3 = 273.16  # [K]      Triple point of water at 611hPa
    # Auxiliary constants
    rdv = 0.6219583540985028  # [ ] rd / rv
    vtmpc1 = 0.6078246934225193  # [ ] rv / rd - 1.0
    vtmpc2 = 0.8608257684344642  # [ ] cpv / cpd - 1.0
    rcpv = -0.46260417446749336  # [ ] cpd / cpv - 1.0
    alvdcp = 2489.2498805542286  # [K] alv / cpd
    alsdcp = 2821.408663799968  # [K] als / cpd
    rcpd = 0.000995381430164039  # [K*kg/J] 1.0 / cpd
    rcvd = 0.0013935340022296545  # [K*kg/J] 1.0 / cvd
    rcpl = 3.1733  # cp_d / cp_l - 1
    clw = 4192.6641119999995  # specific heat capacity of liquid water (rcpl + 1.0) * cpd
    cv_v = 78.37934216297742  # (rcpv + 1.0) * cpd - rv


class IndexConsts(ta.wpfloat, enum.Enum):
    prefactor_r = 14.58
    exponent_r = 0.111
    offset_r = 1.0e-12
    prefactor_i = 1.25
    exponent_i = 0.160
    offset_i = 1.0e-12
    prefactor_s = 57.80
    exponent_s = 0.16666666666666666
    offset_s = 1.0e-12
    prefactor_g = 12.24
    exponent_g = 0.217
    offset_g = 1.0e-08
```

Pattern: subclass `(ta.wpfloat, enum.Enum)` — each member IS a wpfloat (usable directly in arithmetic inside field_operators), grouped by physics domain, precomputed derived values stored as literals with derivation in trailing comments.

---

## 2. Frozen-dataclass config pattern — `SingleMomentSixClassIconGraupelConfig`

File: `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/model/atmosphere/subgrid_scale_physics/microphysics/src/icon4py/model/atmosphere/subgrid_scale_physics/microphysics/single_moment_six_class_gscp_graupel.py`, lines 41–109 (includes `from_fortran_dict`):

```python
@dataclasses.dataclass(frozen=True)
class SingleMomentSixClassIconGraupelConfig:
    """
    Contains necessary parameters to configure icon graupel microphysics scheme.

    Encapsulates namelist parameters.
    Values should be read from configuration.
    Default values are taken from the defaults in the corresponding ICON Fortran namelist files.

    lsedi_ice (option for whether ice sedimendation is performed), lstickeff (option for different empirical formulae of sticking efficiency), and lred_depgrow (option for reduced depositional growth) are removed because they are set to True in the original code gscp_graupel.f90.

    lpres_pri is removed because it is a duplicated option for lsedi_ice.

    ldiag_ttend, and ldiag_qtend are removed because default outputs in icon4py physics granules include tendencies.

    l_cv is removed. This option is to swtich between whether the microphysical processes are isochoric or isobaric. Originally defined as  as input to the graupel in ICON. It is hardcoded to True in ICON.
    The reason for its existence and isobaric may have a positive effect on reducing sound waves. In case it needs to be restored, is_isochoric is a better name.

    The COSMO microphysics documentation "A Description of the Nonhydrostatic Regional COSMO-Model Part II Physical Parameterizations" can be downloaded via this link: https://www.cosmo-model.org/content/model/cosmo/coreDocumentation/cosmo_physics_6.00.pdf
    """

    #: liquid auto conversion mode. Originally defined as isnow_n0temp (PARAMETER) in gscp_data.f90 in ICON. I keep it because I think the choice depends on resolution.
    liquid_autoconversion_option: mphys_options.LiquidAutoConversionType = (
        mphys_options.LiquidAutoConversionType.SEIFERT_BEHENG
    )
    #: snow size distribution interception parameter. Originally defined as isnow_n0temp (PARAMETER) in gscp_data.f90 in ICON. I keep it because I think the choice depends on resolution.
    snow_intercept_option: mphys_options.SnowInterceptParameterization = (
        mphys_options.SnowInterceptParameterization.FIELD_GENERAL_MOMENT_ESTIMATION
    )
    #: Do latent heat nudging. Originally defined as dass_lhn in mo_run_config.f90 in ICON.
    do_latent_heat_nudging: bool = False
    #: Whether a fixed latent heat capacities are used for water. Originally defined as ithermo_water in mo_nwp_tuning_config.f90 in ICON (0 means True).
    use_constant_latent_heat: bool = True
    #: First parameter in RHS of eq. 5.163 in the COSMO microphysics documentation for the sticking efficiency when lstickeff = True (repricated in icon4py because it is always True in ICON). Originally defined as tune_zceff_min in mo_tuning_nwp_config.f90 in ICON.
    ice_stickeff_min: ta.wpfloat = 0.075
    #: Power law coefficient in v-qi ice terminal velocity-mixing ratio relationship, see eq. 5.169 in the COSMO microphysics documentation. Originally defined as tune_zvz0i in mo_tuning_nwp_config.f90 in ICON.
    power_law_coeff_for_ice_mean_fall_speed: ta.wpfloat = 1.25
    #: Exponent of the density factor in ice terminal velocity equation to account for density (air thermodynamic state) change. Originally defined as tune_icesedi_exp in mo_tuning_nwp_config.f90 in ICON.
    exponent_for_density_factor_in_ice_sedimentation: ta.wpfloat = 0.33
    #: Power law coefficient in v-D snow terminal velocity-Diameter relationship, see eqs. 5.57 (depricated after COSMO 3.0) and unnumbered eq. (v = 25 D^0.5) below eq. 5.159 in the COSMO microphysics documentation. Originally defined as tune_v0snow in mo_tuning_nwp_config.f90 in ICON.
    power_law_coeff_for_snow_fall_speed: ta.wpfloat = 20.0
    #: mu exponential factor in gamma distribution of rain particles. Originally defined as mu_rain in mo_nwp_tuning_config.f90 in ICON.
    rain_mu: ta.wpfloat = 0.0
    #: Interception parameter in gamma distribution of rain particles. Originally defined as rain_n0_factor in mo_nwp_tuning_config.f90 in ICON.
    rain_n0: ta.wpfloat = 1.0
    #: coefficient for snow-graupel conversion by riming. Originally defined as csg in mo_nwp_tuning_config.f90 in ICON.
    snow2graupel_riming_coeff: ta.wpfloat = 0.5

    @classmethod
    def from_fortran_dict(
        cls, atmo_dict: dict[str, Any], **overrides: Any
    ) -> SingleMomentSixClassIconGraupelConfig:
        run_nml = atmo_dict["run_nml"]

        nwp_phy_nml = atmo_dict["nwp_phy_nml"]
        nwp_tuning_nml = atmo_dict["nwp_tuning_nml"]
        return cls(
            do_latent_heat_nudging=run_nml["ldass_lhn"],
            use_constant_latent_heat=fortran_config.list_to_value(nwp_phy_nml["ithermo_water"])
            == 0,
            ice_stickeff_min=nwp_tuning_nml["tune_zceff_min"],
            power_law_coeff_for_ice_mean_fall_speed=nwp_tuning_nml["tune_zvz0i"],
            exponent_for_density_factor_in_ice_sedimentation=nwp_tuning_nml["tune_icesedi_exp"],
            power_law_coeff_for_snow_fall_speed=nwp_tuning_nml["tune_v0snow"],
            rain_mu=nwp_phy_nml["mu_rain"],
            rain_n0=nwp_phy_nml["rain_n0_factor"],
            snow2graupel_riming_coeff=nwp_tuning_nml["tune_zcsg"],
            **overrides,
        )
```

Idioms: `@dataclasses.dataclass(frozen=True)`; `#:` sphinx-style attribute comments citing the ICON Fortran origin of every field; typed with `ta.wpfloat` / enums from an `options` module; `from_fortran_dict(cls, atmo_dict, **overrides)` classmethod mapping namelist keys, with `fortran_config.list_to_value` for list-valued namelist entries. Relevant imports at top of file: `from icon4py.model.common.utils import data_allocation as data_alloc, fortran_config`.

---

## 3. Packaged-data conventions

**No icon4py package ships binary data files.** Verified:
- No `.npz` or `.nc` files anywhere under any `src/` tree in `model/`.
- No `data/` directories under any `src/` tree.
- No `importlib.resources` / `importlib_resources` usage anywhere in `model/` source.
- Every package's `package-data` block is the same boilerplate (amps and muphys are identical, `pyproject.toml` line ~79):

```toml
# -- setuptools --
[tool.setuptools.package-data]
'*' = ['*.in', '*.md', '*.rst', '*.txt', 'LICENSE', 'py.typed']

[tool.setuptools.packages]
find = {namespaces = true, where = ['src']}
```

So shipping `.npz` tables inside `src/` in M1 would be a **first** for the repo; the existing precedent for external data is runtime loading from user-supplied paths, as in the muphys driver.

**muphys driver NetCDF loading** — `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/model/atmosphere/subgrid_scale_physics/muphys/src/icon4py/model/atmosphere/subgrid_scale_physics/muphys/driver/common.py` load functions verbatim:

```python
def _as_field_from_nc(
    dataset: netCDF4.Dataset,
    allocator: gtx_typing.Allocator,
    varname: str,
    optional: bool = False,
    dtype: np.dtype | None = None,
) -> gtx.Field[dims.CellDim, dims.KDim] | None:
    if optional and varname not in dataset.variables:
        return None

    var = dataset.variables[varname]
    if var.dimensions[0] == "time":
        var = var[0, :, :]
    data = np.transpose(var)
    if dtype is not None:
        data = data.astype(dtype)
    return gtx.as_field(
        (dims.CellDim, dims.KDim),
        data,
        allocator=allocator,
    )
```

```python
    @classmethod
    def load(
        cls,
        filename: pathlib.Path | str,
        allocator: gtx_typing.Allocator,
        dtype=np.float64,
    ) -> GraupelInput:
        with netCDF4.Dataset(filename, mode="r") as ncfile:
            try:
                ncells = len(ncfile.dimensions["cell"])
            except KeyError:
                ncells = len(ncfile.dimensions["ncells"])

            nlev = len(ncfile.dimensions["height"])

            dz = _calc_dz(np.asarray(ncfile.variables["zg"]).astype(dtype))

            field_from_nc = functools.partial(_as_field_from_nc, ncfile, allocator, dtype=dtype)
            return cls(
                ncells=ncells,
                nlev=nlev,
                dz=gtx.as_field(
                    (dims.CellDim, dims.KDim), np.transpose(dz), allocator=allocator, dtype=dtype
                ),
                t=field_from_nc("ta"),
                p=field_from_nc("pfull"),
                qs=field_from_nc("qs"),
                qi=field_from_nc("cli"),
                qg=field_from_nc("qg"),
                qv=field_from_nc("hus"),
                qc=field_from_nc("clw"),
                qr=field_from_nc("qr"),
                rho=field_from_nc("rho"),
            )
```

(`GraupelOutput.load` follows the same shape with `optional=True` for precipitation fields; `GraupelOutput.allocate` is the field-allocation counterpart using `functools.partial(gtx.zeros, domain=..., allocator=..., dtype=...)`.) Pattern: dataclass container + `load(filename, allocator, dtype)` classmethod + `functools.partial` over a private `_as_field_from_nc` helper; `filename` is always a caller-supplied path, never a packaged resource.

---

## 4. `StencilTest` base class

File: `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/model/testing/src/icon4py/model/testing/stencil_tests.py`. Supporting module-level `Output` dataclass:

```python
@dataclasses.dataclass(frozen=True)
class Output:
    name: str
    refslice: tuple[slice, ...] = dataclasses.field(default_factory=lambda: (slice(None),))
    gtslice: tuple[slice, ...] = dataclasses.field(default_factory=lambda: (slice(None),))
```

The class in full (lines 179–313):

```python
class StencilTest:
    """
    Base class to be used for testing stencils.

    Example (pseudo-code):

        >>> class TestMultiplyByTwo(StencilTest):  # doctest: +SKIP
        ...     PROGRAM = multiply_by_two  # noqa: F821
        ...     OUTPUTS = ("some_output",)
        ...     STATIC_PARAMS = {"category_a": ["flag0"], "category_b": ["flag0", "flag1"]}
        ...
        ...     @pytest.fixture
        ...     def input_data(self):
        ...         return {"some_input": ..., "some_output": ...}
        ...
        ...     @staticmethod
        ...     def reference(some_input, **kwargs):
        ...         return dict(some_output=np.asarray(some_input) * 2)
    """

    PROGRAM: ClassVar[gtx_typing.Program | gtx_typing.FieldOperator]
    OUTPUTS: ClassVar[tuple[str | Output, ...]]
    STATIC_PARAMS: ClassVar[dict[str, Sequence[str]] | None] = None

    reference: ClassVar[Callable[..., Mapping[str, np.ndarray | tuple[np.ndarray, ...]]]]

    @pytest.fixture
    def _configured_program(
        self,
        backend_like: model_backends.BackendLike,
        static_variant: Sequence[str],
        input_data: dict[str, gtx.Field | tuple[gtx.Field, ...]],
        grid: base.Grid,
    ) -> Callable[..., None]:
        unused_static_params = set(static_variant) - set(input_data.keys())
        if unused_static_params:
            raise ValueError(
                f"Parameter defined in 'STATIC_PARAMS' not in 'input_data': {unused_static_params}"
            )
        static_args = {name: [input_data[name]] for name in static_variant}
        backend = model_options.customize_backend(self.PROGRAM, backend_like)
        program = self.PROGRAM.with_backend(backend)
        if backend is not None:
            if isinstance(program, FieldOperator):
                if len(static_args) > 0:
                    raise NotImplementedError(
                        "'FieldOperator's do not support static arguments yet."
                    )
            else:
                program.compile(
                    offset_provider=grid.connectivities,
                    **static_args,  # type: ignore[arg-type]
                )

        test_func = device_utils.synchronized_function(program, allocator=backend)
        return test_func

    @pytest.fixture
    def _properly_allocated_input_data(
        self,
        input_data: dict[str, gtx.Field | tuple[gtx.Field, ...]],
        backend_like: model_backends.BackendLike,
    ) -> dict[str, Any]:
        # TODO(havogt): this is a workaround,
        # because in the `input_data` fixture provided by the user
        # it does not allocate for the correct device.
        allocator = model_backends.get_allocator(backend_like)
        return allocate_data(allocator=allocator, input_data=input_data)

    def _verify_stencil_test(
        self,
        input_data: dict[str, gtx.Field | tuple[gtx.Field, ...]],
        reference_outputs: Mapping[str, np.ndarray | tuple[np.ndarray, ...]],
    ) -> None:
        for out in self.OUTPUTS:
            name, refslice, gtslice = (
                (out.name, out.refslice, out.gtslice)
                if isinstance(out, Output)
                else (out, (slice(None),), (slice(None),))
            )

            input_data_name = input_data[name]  # for mypy
            # TODO(iomaganaris, havogt, nfarabullini): tolerance was increased from 1e-7 to 1e-6
            # to cover floating point descripancies observed in CI tests. Failing CI can be found in
            # https://gitlab.com/cscs-ci/ci-testing/webhook-ci/mirrors/5125340235196978/2255149825504673/-/pipelines/2184694383
            # from PR#861. Reason is probably derivatives of random data. Investigate and lower tolerance back to 1e-7 if possible.
            relative_tolerance = 3e-6
            if isinstance(input_data_name, tuple):
                for i_out_field, out_field in enumerate(input_data_name):
                    test_utils.assert_dallclose(
                        out_field.asnumpy()[gtslice],
                        reference_outputs[name][i_out_field][refslice],
                        equal_nan=True,
                        err_msg=f"Verification failed for '{name}[{i_out_field}]'",
                        rtol=relative_tolerance,  # TODO(iomaganaris, havogt, nfarabullini): check above comment
                    )
            else:
                reference_outputs_name = reference_outputs[name]  # for mypy
                assert isinstance(reference_outputs_name, np.ndarray)
                test_utils.assert_dallclose(
                    input_data_name.asnumpy()[gtslice],
                    reference_outputs_name[refslice],
                    equal_nan=True,
                    err_msg=f"Verification failed for '{name}'",
                    rtol=relative_tolerance,  # TODO(iomaganaris, havogt, nfarabullini): check above comment
                )

    @staticmethod
    def static_variant(request: pytest.FixtureRequest) -> Sequence[str]:
        """
        Fixture for parametrization over the `STATIC_PARAMS` of the test class.

        Note: the actual `pytest.fixture()`  decoration happens inside `__init_subclass__`,
          when all information is available.
        """
        _, variant = request.param
        return () if variant is None else variant

    def __init_subclass__(cls, **kwargs: Any) -> None:
        super().__init_subclass__(**kwargs)

        setattr(cls, f"test_{cls.__name__}", test_and_benchmark)

        # decorate `static_variant` with parametrized fixtures, since the
        # parametrization is only available in the concrete subclass definition
        if cls.STATIC_PARAMS is None:
            # not parametrized, return an empty tuple
            cls.static_variant = staticmethod(pytest.fixture(lambda: ()))  # type: ignore[method-assign] # we override with a non-parametrized function
        else:
            cls.static_variant = staticmethod(  # type: ignore[method-assign]
                pytest.fixture(params=cls.STATIC_PARAMS.items(), scope="class", ids=lambda p: p[0])(
                    cls.static_variant
                )
            )
```

Key mechanics: subclass sets `PROGRAM`, `OUTPUTS`, a `reference` staticmethod, and an `input_data` pytest fixture; `__init_subclass__` injects the actual test function (`test_and_benchmark(self, benchmark, grid, _properly_allocated_input_data, _configured_program, request)`, module-level, lines 84–176) as `test_<ClassName>`. The injected test requires `grid` and `backend_like` fixtures to exist — supplied via conftest (see 5). Verification tolerance is fixed at `rtol=3e-6` through `test_utils.assert_dallclose`.

---

## 5. pytest markers/fixtures for tests

**muphys stencil_tests conftest** (`/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/model/atmosphere/subgrid_scale_physics/muphys/tests/muphys/stencil_tests/conftest.py`) — the entire wiring is two imports:

```python
from icon4py.model.testing.fixtures.datatest import backend_like
from icon4py.model.testing.fixtures.stencil_tests import grid, grid_manager
```

(The package-level `tests/muphys/conftest.py` is empty apart from the license header. `model/common/tests/common/math/stencil_tests/conftest.py` is identical two-line wiring.) The `grid` fixture defaults to `DEFAULT_GRID = "simple"` with `DEFAULT_NUM_LEVELS = 40` (from `model/testing/src/icon4py/model/testing/fixtures/stencil_tests.py`), so StencilTest subclasses run out-of-the-box without data downloads.

**Root pytest config** (`/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/pyproject.toml`):

```toml
[tool.pytest.ini_options]
addopts = ['-p icon4py.model.testing.pytest_hooks', '--strict-markers']
markers = [
  "embedded_remap_error",
  "embedded_static_args",
  "skip_value_error",
  "datatest",
  "embedded_only",
  "uses_concat_where",
  "cpu_only",
  "infinite_concat_where",
  "gtfn_too_slow",
  "continuous_benchmarking",
  "benchmark_only: benchmark only tests without verification"
]
```

`--strict-markers` means any new marker must be registered here. `pythonpath` already includes `model/atmosphere/subgrid_scale_physics/amps`.

**Pure-numpy unit tests need NO fixtures.** Precedent: `model/common/tests/common/io/unit_tests/test_cf_utils.py` has no conftest anywhere on its path (no conftest exists in `model/common/tests/common/io/` at all) and is plain pytest:

```python
@pytest.mark.parametrize("input_", test_io_utils.state_values())
def test_to_canonical_dim_order(input_):
    input_dims = input_.dims
    output = cf_utils.to_canonical_dim_order(input_)
    assert output.dims == (input_dims[1], input_dims[0])
```

And the amps package already has one: `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/model/atmosphere/subgrid_scale_physics/amps/tests/amps/unit_tests/test_package.py`:

```python
def test_package_imports_and_has_version():
    from icon4py.model.atmosphere.subgrid_scale_physics import amps  # noqa: PLC0415

    assert amps.__version__ == "0.2.0"
```

Conclusion: plain pytest is fine for non-grid unit tests; the `grid`/`backend_like` conftest imports are needed only in directories containing `StencilTest` subclasses.

---

## 6. amps spikes — docstrings + adopted-idiom constants

All in `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/model/atmosphere/subgrid_scale_physics/amps/spikes/`.

### spikes/common.py (shared helpers; module docstring + constants)

```python
"""Shared helpers for M0 feasibility spikes. Not part of the amps package API."""
```
```python
NCELLS = 4096
NLEV = 61
GENDIR = pathlib.Path(__file__).parent / "_generated"
```
Also codifies: `backends()` returns `{"embedded": None, "gtfn_cpu": run_gtfn_cached}` (from `gt4py.next.program_processors.runners.gtfn import run_gtfn_cached`), and `load_generated_operator(source, module_name, attr)` writes generated DSL source to a real file "(gt4py parses inspect.getsource)" before importing it.

### spike_a_remap_gather.py

```python
"""Spike A: gather-formulated bin remap and value-indexed table lookup via
as_offset. Answers: does as_offset(Koff, int_field) work on embedded and
gtfn_cpu, for (a) field self-gather with a computed shift, (b) a K-stored
table gathered at an index computed from field values?

Follow-up (same task): the plain K-only-table idiom below is a NO-GO (see the
comment block under `_table_gather`, kept for the record). Two rescue variants
were probed: broadcasting the table to (Cell, K) before the gather (still a
NO-GO, differently on each backend) and passing the table pre-tiled to full
(Cell, K) shape, i.e. memory-replicated across cells (a GO on both backends).

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_a_remap_gather.py
"""
```
No module-level constants. Adopted idiom: **tiled (Cell, K) tables** (`_table_gather_tiled`) — the only both-backend GO for value-indexed table lookup.

### spike_b_collection_codegen.py

```python
"""Spike B: code-generate an unrolled bin-pair collection (coalescence) kernel
and measure embedded/gtfn_cpu compile + run cost at nbins = 8, 20, 40.

Physics stand-in: Golovin kernel K_ij = k0*(m_i+m_j), mass-doubling bins
m_b = 2^b, Kovetz-Olund style two-bin deposit of coalesced mass, plus loss
terms. Structure (nbins^2 pair terms scattered to <=2 destination bins each,
plus nbins x nbins loss sums) matches the real AMPS coalescence shape.

This is the compile-time gate for the whole AMPS port architecture (spec §9
risk 1): can a single generated, fully-unrolled field_operator for the 40-bin
case even compile on gtfn_cpu in a bounded time? A backend failure or an
excessive compile time, captured verbatim, is itself the measurement.

CRITICAL for measurement validity: the environment has a PERSISTENT on-disk
gt4py build cache (GT4PY_BUILD_CACHE_LIFETIME=PERSISTENT). `clear_gt4py_cache`
below removes the worktree-local `.gt4py_cache` directory (the cache gt4py
resolves to when GT4PY_BUILD_CACHE_DIR is unset and the process cwd is the
worktree root, as used by the documented run command) before each gtfn_cpu
compile so first-call timings are genuine cold-compile numbers, not artifacts
of a warm on-disk cache from a previous run/spike.

WARNING: with the muphys-precedented recursion-limit raise below, nbins=40's
gtfn_cpu compile SUCCEEDS but takes ~2579s (~43 minutes, measured) instead of
crashing in ~13s. Running this script end-to-end via the documented command
will therefore take on the order of 45-50 minutes total (dominated by that
one compile). This is an intentional, verified-necessary trade: without the
raised limit, nbins=40 crashes with RecursionError instead (see the task-8
report for both numbers). There is no way to get nbins=40's real gtfn_cpu
number without paying one or the other.

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_b_collection_codegen.py
"""
```
Module constants:
```python
RTOL_BY_NBINS: dict[int, float] = {8: 1e-12, 20: 1e-12, 40: 1e-6}


HEADER = """\
import gt4py.next as gtx

from icon4py.model.common import field_type_aliases as fa, type_alias as ta

"""
```
Recursion-limit idiom (line 232): `with muphys_driver_utils.recursion_limit(10**5):` where the import is `from icon4py.model.atmosphere.subgrid_scale_physics.muphys.driver import utils as muphys_driver_utils`.

### spike_c_wide_scan.py

```python
"""Spike C: scan_operator with an 80-field NamedTuple carry (2x40, see "FIX"
below) -- the shape of 40-bin implicit sedimentation (muphys does 4
hydrometeor classes; AMPS needs 40 bins). Generated source (explicit
80-field NamedTuple + 40 scalar args).

Numerics stand-in (per level, forward/downward scan, muphys-style):
  flux_b_new = v_b * q_b * rho ; q_out_b = q_b + zeta*(carry.f_b - flux_b_new)
carry = (fluxes, updated q). Correctness checked vs numpy recurrence.
Perturbation check (run_perturbation_check, always run): a level-0
perturbation must propagate into level 1's output, proving this is a
genuine sequential scan and not a degenerate level-parallel computation.

FIX (post-review): the first version of this spike (commit c2702ec60) had a
degenerate carry -- the exposed output (`fl_b`) never read `carry` at all
(only the *unused* `qo_b` intermediate did), so the whole "scan" was
mathematically 40 independent per-level pointwise computations with no
actual cross-level data dependency, and its 82.7s gtfn_cpu compile number
measured that degenerate case, not a genuine wide-carry scan. Fixed by:
widening the carry to 80 fields (f_b: flux, q_b: updated q -- q_b is
carried but deliberately not read back, per review: this stress-tests
carry *width* without adding a second live recurrence) and exposing
`s.q_XX` (which genuinely reads `carry.f_b`) as the field_operator's output
instead of `s.f_XX`. See "MEASURED" below for the corrected number.

CRITICAL for measurement validity: the environment has a PERSISTENT on-disk
gt4py build cache (GT4PY_BUILD_CACHE_LIFETIME=PERSISTENT). `clear_gt4py_cache`
below removes the worktree-local `.gt4py_cache` directory (the cache gt4py
resolves to when GT4PY_BUILD_CACHE_DIR is unset and the process cwd is the
worktree root, as used by the documented run command) before the gtfn_cpu
compile so the first-call timing is a genuine cold-compile number, not an
artifact of a warm on-disk cache from a previous run/spike (same idiom as
spike_b_collection_codegen.py; the helper functions are copied verbatim
below since spikes are standalone scripts).

SCAN OUTPUT CONVENTION (measured, not assumed): gt4py 1.1.11's scan_operator
emits the *post-update* carry at each level -- the NamedTuple a scan step
*returns* is both next level's input carry AND the value written to this
level's output field. numpy_reference below implements this convention
(prev_flux carried in, this level's fl computed and stashed for next level,
this level's *output* is the carry-dependent qo). See run_perturbation_check
for the executable proof this is genuinely sequential.

MEASURED (this spike's actual go/no-go datum, corrected/superseded number --
see FIX above; the original 82.7s measured a degenerate, non-sequential
computation): gtfn_cpu -- clear-cache cold compile succeeds in 186.9s,
steady-state run 115.9ms/call (NBINS=40 i.e. an 80-field carry, NCELLS=4096,
NLEV=61; recursion_limit(10**5) applied per the muphys driver precedent).
Still a clean GO by spike_b_collection_codegen.py's thresholds (<120s clean
go, 120s-15min go-with-caching, >15min/crash no-go) -- 186.9s lands just
past the "clean" boundary into "go with per-kernel compile caching" (which
run_gtfn_cached already provides within a process), and is ~2.3x the
degenerate 82.7s number, consistent with the genuinely sequential carry
(now threading `carry.f_b` into the exposed output, plus double the carry
field count) costing more to compile/lower than the old level-parallel
computation, while still being ~14x faster than spike_b's O(nbins^2)
collection kernel at the same nbins=40 (2578.7s).

embedded is NOT attempted at full scale by default: a prior full-scale
(NCELLS=4096/NLEV=61) run of this spike's (pre-fix) generated operator ran
for ~90-95 minutes at ~100% CPU and then the process was simply gone -- no
RESULT line, no Python traceback (stderr was merged into the same captured
stream), nothing; not confirmed via OS logs, but consistent with an
out-of-memory kill. See the task-9 report for the full timeline. Given that
observed cost/risk, embedded now defaults to a small-NCELLS calibration
ladder (reproducible from this committed script, see EMBEDDED_CALIBRATION_NCELLS)
instead of the full documented grid size; set AMPS_SPIKE_C_FULL_EMBEDDED=1 to
attempt the full common.NCELLS scale anyway (not recommended without a
process-level time/memory guard).

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_c_wide_scan.py
(gtfn_cpu now runs first, by default, followed by the small-NCELLS embedded
calibration ladder. AMPS_SPIKE_C_GTFN_ONLY=1 / AMPS_SPIKE_C_EMBEDDED_ONLY=1
isolate one backend; AMPS_SPIKE_C_FULL_EMBEDDED=1 opts into the
(discouraged, see above) full-scale embedded attempt.)
"""
```
Module constants:
```python
NBINS = 40
VTS = np.linspace(0.1, 8.0, NBINS)  # per-bin fall speeds, folded as constants

# embedded default: a small-NCELLS calibration ladder (see module docstring
# for why full-scale embedded is opt-in, not default). Chosen to match the
# range independently calibrated during the original (pre-fix) spike run
# (linear ~55-56 ms/cell at these sizes for the degenerate carry; re-checked
# below for the fixed, genuinely-sequential carry).
EMBEDDED_CALIBRATION_NCELLS = (16, 64, 256)
```
Recursion-limit idiom (line 311): `with muphys_driver_utils.recursion_limit(10**5):`.

### spike_d_esat.py

```python
"""Spike D: saturation vapor pressure over liquid -- analytic Murphy-Koop in
DSL vs Fortran-style table + linear interpolation (estbar: 150 entries,
T index = int(T)-163, clamped). Decides whether the port replaces the
Fortran LUTs with analytic formulas (spec 4, core/thermo.py).

Fortran reference: estbar(i) tabulates Murphy-Koop at T = 163+i K
(i = 1..150), linear interp between entries; see QSPARM2 in
mod_amps_utility.F90 (scale_amps repo).

Relationship to Spike A (spike_a_remap_gather.py): that spike found the
*plain* K-only-table-gather idiom (`table(as_offset(Koff, expr))` returned
directly, with no further arithmetic) is a decoration-time NO-GO -- FOAST
return-type deduction derives the result dims purely from the remapped
field's own dims (K-only for a KDim-only table), ignoring the extra CellDim
carried by the *offset expression*, so the annotated `CellKField` return
type never matches and `@gtx.field_operator` raises `DSLError` before any
backend runs.

The brief for this spike assumed `_esat_table_k_only` below (structurally:
K-only table, CellK-computed offset) would hit that exact same decoration
failure. It does NOT -- verified by actually attempting it (not assumed from
Spike A's finding). The two cases differ in one relevant way: Spike A's
`_table_gather` *returns the gathered field directly*; `_esat_table_k_only`
instead combines *two* gathered fields (`e0`, `e1`) arithmetically with
`frac` (a genuinely `(Cell, K)`-shaped intermediate, since `frac` derives
from `t`). That extra arithmetic changes which code path FOAST's dims
deduction takes:

  - Decoration: SUCCEEDS (contra the brief's expectation) -- the annotated
    `CellKField` return type is accepted.
  - embedded execution: FAILS at run time (not decoration time), for both
    `offset_provider={}` and `offset_provider={"Koff": dims.KDim}`, with:
      ValueError("Dimensions 'K[vertical], Cell[horizontal]' are not
      ordered correctly, expected 'Cell[horizontal], K[vertical]'.")
    i.e. the deduced dims *do* include both Cell and K (unlike Spike A's
    pure-K-only deduction), but in the wrong order, and gt4py's embedded
    executor enforces canonical (Cell, K) ordering at run time rather than
    at decoration time.
  - gtfn_cpu execution: unexpectedly SUCCEEDS with `offset_provider={}` and
    produces numerically plausible values (max relative error vs analytic
    ~4e-3 at a 16-cell probe scale, consistent with 1 K table resolution) --
    see `run_k_only_variant` below for the full-scale, reproducible number.
    This asymmetry (gtfn_cpu tolerating a dims-order mismatch that embedded
    rejects) was not observable in Spike A, since that spike's K-only case
    never got past decoration on either backend.

Net: the K-only table idiom is *not* a clean decoration-time NO-GO here, but
it is backend-inconsistent (fails on embedded, "works" on gtfn_cpu) and thus
not adopted. `_esat_table_tiled` below -- the memory-replicated fallback
Spike A proved is a GO on both backends -- is the measured, portable path;
its RESULT lines gate this spike's correctness asserts and timings.

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_d_esat.py
"""
```
Module constant:
```python
NLEV = 200  # table needs >= 150 K-levels for the gather idiom
```

### spike_e_counter_rng.py

```python
"""Spike E: counter-based RNG expressible in the GT4Py DSL, which allows only
+ - * // % on integers (no bitwise ops, no shifts -- verified against gt4py
1.1.11 func_to_foast). Construction: an int64 counter mixed from
(cell, k, bin, step), one nonlinear (quadratic) mixing round, then 3
Park-Miller Lehmer rounds modulo the Mersenne prime 2^31-1.

Quality bar (for habit selection, NOT crypto): mean ~ 0.5, var ~ 1/12,
lag-1 correlation across each axis < 0.01, coarse 16-bin histogram flat
within 1%. All checked vs numpy replica (bit-identical integers) -- that
match is the gating assert for everything downstream.

Two findings diverge from the brief; both driven by actually running things,
not assumed (full numeric trail in task-11-report.md):

FINDING 1 (statistical construction): the brief's 3-round-Lehmer-only
recipe -- x = combine(cell,k,bin,step) % M31; x = (x*A1) % M31;
x = (x*A2) % M31; x = (x*A3) % M31 -- FAILS the lag-1 correlation bar badly
(lag1_cell=-0.337, lag1_k=-0.496 vs the <0.01 bar). The brief's prescribed
remedy for a failing assert ("add one more LCG round, A4=40692") does NOT
fix it (lag1_cell=-0.387, lag1_k=-0.442 with 4 rounds); rounds 5-7 (tried
out of caution) also fail, with lag1 in [0.27, 0.82] and no trend toward
zero. Root cause: `combine` is affine in (cell, k, bin, step) mod M31, and
each Lehmer round x -> (x*A) % M31 is a linear map over the ring Z_M31;
composing linear maps with an affine one stays affine. Along a unit-step
axis (consecutive cell or k indices), an affine map mod M31 is exactly a
Weyl/sawtooth sequence x_n = n*delta + c (mod M31) -- its lag-1 correlation
is a deterministic function of the single scalar slope `delta` and is not
damped by stacking more multiplicative rounds, since composing more linear
maps only changes `delta`, never removes the affine structure that causes
the correlation. This is structural, not statistical noise curable by "one
more round" of the same kind. FIX: insert one genuinely nonlinear mixing
round between `combine` and the Lehmer rounds -- squaring the running
state, `x = (x*x + x + 1) % M31`. Squaring makes the state quadratic (not
affine) in the original index, breaking the Weyl structure; `x*x`, `+`, `%`
are all DSL-legal integer ops, no new primitives needed. `x` stays < M31 ~=
2.1e9 at every step, so `x*x` peaks at ~4.6e18 -- comfortably inside
int64's ~9.22e18 max, no overflow. This one change takes all four quality
metrics from a hard fail to a clean pass (see RESULT rng_quality line),
still with only the brief's original 3 Lehmer rounds -- final construction
is 1 quadratic round + 3 Lehmer rounds (4 total mixing rounds after the
initial affine combine); no A4 needed.

FINDING 2 (DSL API, the potential pitfall the task brief flagged): the
brief's `_hash01` references its mixing multipliers (C_CELL, M31, A1, ...)
as bare module-level Python int globals inside the field_operator body.
That fails two different ways, verified in isolation before touching this
file:
  (a) a bare Python int literal in an int64 expression defaults to int32
      and gt4py rejects the implicit widen: `DSLError: Could not promote
      'Field[[Cell], int64]' and 'int32' to common type in call to '*'`.
  (b) wrapping the module constant as np.int64(...) fixes (a) (embedded
      then runs), but gtfn_cpu's ITIR lowering fails to resolve the
      constant as a closure symbol at all: `EveValueError: Symbols
      {SymbolRef('M31'), SymbolRef('C_CELL'), ...} not found`. This
      reproduces regardless of whether the reference is bare or wrapped in
      `astype(NAME, gtx.int64)` -- the closure *variable name* itself is
      what gtfn_cpu fails to resolve, not the type promotion.
  A bare integer *literal* (not a named global) written directly in the
  operator body, e.g. `astype(999983, gtx.int64)`, does resolve correctly
  on both backends -- but that would mean hardcoding all eight mixing
  constants as unnamed magic numbers in the body, which is worse than
  either failure mode above. The adopted fix instead promotes every mixing
  constant (c_cell, c_k, c_bin, c_step, m31, a1, a2, a3) to an explicit
  `gtx.int64` scalar parameter of `_hash01`, passed at each call site as
  `gtx.int64(C_CELL)` etc. from the still-named module-level Python int
  constants -- scalar call arguments are threaded through the compiled
  program directly and hit neither failure mode; `bin_id`/`step` were
  already scalar params in the brief, so this generalizes that existing
  idiom rather than introducing a new one. `gtx.int64(7)` scalar
  construction at the call site (the brief's other flagged pitfall) worked
  with no issue on either backend.

Run: uv run --frozen python model/atmosphere/subgrid_scale_physics/amps/spikes/spike_e_counter_rng.py
"""
```
RNG construction constants (module-level):
```python
M31 = 2147483647  # 2^31 - 1 (Mersenne prime, Park-Miller modulus)
A1 = 16807  # Park-Miller multipliers for the (linear) Lehmer rounds
A2 = 48271
A3 = 69621
# mixing multipliers for the initial counter combine (distinct odd primes,
# < 2^20 to keep products of int32-range counters safely inside int64)
C_CELL = 999983
C_K = 424243
C_BIN = 786433
C_STEP = 611953
```

**recursion_limit precedent** — `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/model/atmosphere/subgrid_scale_physics/muphys/src/icon4py/model/atmosphere/subgrid_scale_physics/muphys/driver/utils.py`:

```python
@contextlib.contextmanager
def recursion_limit(limit: int) -> Generator[None, None, None]:
    original_limit = sys.getrecursionlimit()
    sys.setrecursionlimit(limit)
    try:
        yield
    finally:
        sys.setrecursionlimit(original_limit)
```

muphys's own drivers use `utils.recursion_limit(10**4)` (`run_graupel_only.py:80`, with `# TODO(havogt): make an option in gt4py?`) and `10**5` (`run_full_muphys.py:121,149`); spikes b and c import it as `muphys_driver_utils` and use `10**5`.

---

## 7. tach.toml `[external]` exclude block

File: `/Users/jcanton/projects/icon4py/.worktrees/amps_microphysics/tach.toml`, verbatim including the comments:

```toml
# exclude optional external dependencies
[external]
exclude = [
    # workspace member: tach resolves this import via source_roots, so the external pip dep appears unused
    "icon4py_tools",
    # icon4py-atmosphere-amps (scaffold, Task 6) declares these ahead of the GT4Py spikes
    # (Tasks 7-11) that will import them; no source file uses them yet. Remove once they do.
    "icon4py_common",
    "numpy",
    "gt4py",
    "fprettify",
    "configargparse",
    "cupy",
    "ghex",
    "dace",
    "mpi4py",
    "netcdf4",
    "pymetis",
    "xarray",
    "uxarray",
    "cftime",
    "viztracer",
    "clang_format",
]
rename = ["serialbox:serialbox4py"]
```

Per the comment, M1 must remove `icon4py_common`, `numpy`, `gt4py` (etc.) from this list as real imports land in `amps/src/`. Note the amps module dependency declaration also exists (`[[modules]] path = "icon4py.model.atmosphere.subgrid_scale_physics.amps"` around line 51; muphys's declares `depends_on = [{ path = "icon4py.model.common" }]`).