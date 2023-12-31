// pybind11
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

// pydisort
#include <cppdisort/cppdisort.hpp>

namespace py = pybind11;

PYBIND11_MODULE(pydisort, m) {
  m.doc() = R"(
  Python bindings for DISORT (Discrete Ordinates Radiative Transfer) Program.

  Summary
  -------
  This module provides a python interface to the C version of the DISORT program.
  Please consult the DISORT publication [1]_ for more information on the DISORT program,
  and the C-DISORT C publication [2]_ for more information on the C version of the DISORT program.

  Small changes have been made to the C-DISORT program to make it compatible with python scripting.
  The C-DISORT program has been wrapped first in a C++ class (DisortWrapper), and the C++ class
  has been bound to python using pybind11. 

  Pydisort features the following benefits over the original C-DISORT program:

  - Proper handling of errors rather than abrupt exit of the program. Errors
    can be caught and and handled in the python script.
  - Memory management is handled by the C++ class. The user does not need to
    worry about memory allocation and deallocation.
  - Documentation is automated using sphinx and readthedocs.
  - Safety guards are implemented to prevent the user from setting incorrect
    values for arrays or calling methods in the wrong order.

  Note that the underlying calculation engine is still the same as the C-DISORT program. 
  So the speed of pydisort is the same as the origin C-DISORT program.

  Examples
  --------
  - Example 1: Calculate attenuation of radiative flux in a plane-parallel atmosphere

  >>> from pydisort import disort, RFLDIR, FLDN, FLUP
  >>> ds = disort()
  >>> flags = {"onlyfl": True, "usrtau": False, "usrang": False}
  >>> ds.set_flags(flags)
  >>> ds.set_atmosphere_dimension(nlyr = 4)
  >>> ds.seal()
  >>> ds.set_optical_thickness([0.1, 0.2, 0.3, 0.4])
  >>> ds.fbeam = 3.14159
  >>> _, flx = ds.run()
  >>> flx = flx[:, [RFLDIR, FLDN, FLUP]]
  >>> flx
  array([[3.14159   , 0.        , 0.        ],
         [2.84262818, 0.        , 0.        ],
         [2.32734711, 0.        , 0.        ],
         [1.72414115, 0.        , 0.        ],
         [1.15572637, 0.        , 0.        ]])

  In the example above, flx has two dimensions. The first dimension is number of atmosphere
  levels (nlyr + 1 = 5), and the second dimension is number of extracted flux fields (3).
  ``RFDLIR``, ``FLDN``, and ``FLUP`` are the indices representing the three flux fields:
  direct flux, diffuse flux, and upward flux, respectively. 
  Consult ``pydisort.run()`` method for more information on the indices of flux fields.
  The attenuation of radiative flux is acoording to the Beer-Lambert law, i.e.,

  .. math::

    F(z) = F(0) \exp(-\tau(z)),

  where :math:`F(z)` is the radiative flux at level :math:`z`, 
  :math:`F(0)` is the radiative flux at the top of the atmosphere, and :math:`\tau(z)` is the
  optical depth from the top of the atmosphere to level :math:`z`. The default direction of
  radiative flux is nadir.

  - Example 2: Calculate intensity from isotropic scattering in a plane-parallel atmosphere

  >>> from pydisort import disort, get_phase_function
  >>> ds = disort()
  >>> flags = {"planck": False}
  >>> ds.set_flags(flags)
  >>> ds.set_atmosphere_dimension(nlyr = 1, nstr = 16, nmom = 16)
  >>> ds.set_intensity_dimension(nuphi = 1, nutau = 2, numu = 6)
  >>> ds.seal()
  >>> pmom = get_phase_function(nmom = 16, model = "isotropic")
  >>> ds.set_optical_thickness([0.1])
  >>> ds.set_single_scattering_albedo([1.0])
  >>> ds.set_phase_moments(pmom)
  >>> ds.set_user_optical_depth([0.0, 0.1])
  >>> ds.set_user_cosine_polar_angle([-1.0, -0.5, -0.1, 0.1, 0.5, 1.0])
  >>> ds.set_user_azimuthal_angle([0.0])
  >>> ds.umu0 = 1.0
  >>> ds.fbeam = 3.14159
  >>> rad, flx = ds.run()
  >>> rad
  array([[[0.        , 0.        , 0.        , 0.18095504, 0.0516168 , 0.02707849],
          [0.02703935, 0.05146774, 0.17839685, 0.        , 0.        , 0.        ]]])

  The intensity array ``rad`` has three dimensions. The first dimension is the 
  azimuthal angles (1). The second dimension is optical depth (2). 
  The third dimension is the polar angles (6). The result is interpreted as backsattering
  at the top of the atmosphere (optical depth = 0.0) and forward scattering at the bottom
  of the atmosphere (optical depth = 0.1). 
  

  Troubleshooting
  ---------------
  - The most common error is "RuntimeError: DisortWrapper::Run failed". When
    this error occurs, please check the error message printed before the
    error message. The error message printed before the error message
    usually provides more information on the cause of the error. Once you identify
    the cause of the error, you can fix the error by calling ``unseal()`` method,
    then setting the correct values, and then calling ``seal()`` method again.

  - One common issue that results in "RuntimeError: DisortWrapper::Run failed"
    is incompatible flags for flux or intensity calculations. For example, ``usrtau``
    and ``usrang`` flags should set to ``False`` when ``onlyfl`` flag is set to ``True``.

  - Another common issue is setting the wrong values for temperature, optical thickness,
    single scattering albedo, or phase function moments. All these values must be
    positive. If you set a negative value, you will get "RuntimeError: DisortWrapper::Run failed".

  - The program should not exit unexpectedly. If the program exits unexpectedly,
    please report the issue to the author (zoey.zyhu@gmail.com).

  Important Dimensions
  --------------------
  nlyr
      number of atmosphere layers
  nstr
      number of radiation streams
  nmom
      number of phase function moments in addition to the zeroth moment

  Tips
  ----
  - Number of atmosphere levels is one more than the number of atmosphere layers.

  - Temperature is defined on levels, not layers. Other properties such as
    optical thickness, single scattering albedo, and phase function moments
    are defined on layers.

  - You can use ``print()`` method to print some of the DISORT internal states.

  - You can chain methods such as ``set_flags``, ``set_atmosphere_dimension()``, 
    ``set_intensity_dimension()``, and ``seal()`` together.

  - If you want to have more insights into DISORT internal inputs,
    you can set the ``print-input`` flag to ``True``.
    The DISORT internal inputs will be printed to the standard output
    when the ``run()`` method is called.

  References
  ----------
  .. [1] Stamnes, K., Tsay, S. C., Wiscombe, W., & Jayaweera, K. (1988). 
         Numerically stable algorithm for discrete-ordinate-method radiative transfer in multiple scattering and emitting layered media. 
         Applied Optics, 27(12), 2502-2509.
  .. [2] Buras, R., & Dowling, T. (1996). 
         Discrete-ordinate-method for radiative transfer in planetary atmospheres: Generalization of the doubling and adding method. 
         Journal of Quantitative Spectroscopy and Radiative Transfer, 55(6), 761-779.
  )";

  m.attr("RFLDIR") = 0;
  m.attr("FLDN") = 1;
  m.attr("FLUP") = 2;
  m.attr("DFDT") = 3;
  m.attr("UAVG") = 4;
  m.attr("UAVGDN") = 5;
  m.attr("UAVGUP") = 6;
  m.attr("UAVGSO") = 7;

  m.def("get_phase_function", &get_phase_function, R"(
      Get phase function moments based on a phase function model

      Parameters
      ----------
      nmom : int
          Number of phase function moments
      model : str
          Phase function model. 
      gg : float
          Asymmetry factor for henyey-greenstein phase function

      Returns
      -------
      pmom : List[float]
          Phase function moments, shape (1 + nmom,)

      Notes
      -----
      The following phase function models are supported:

      .. list-table::
          :widths: 25 40
          :header-rows: 1

          * - Model
            - Description
          * - 'isotropic'
            - Isotropic phase function, [1, 0, 0, 0, ...]
          * - 'rayleigh'
            - Rayleigh scattering phase function, [1, 0, 0.1, 0, ...]
          * - 'henyey_greenstein'
            - Henyey-Greenstein phase function, [1, gg, gg^2, gg^3, ...]
          * - 'haze_garcia_siewert'
            - Tabulated haze phase function by Garcia/Siewert
          * - 'cloud_garcia_siewert'
            - Tabulated cloud phase function by Garcia/Siewert
      )",
      py::arg("nmom"), py::arg("model"), py::arg("gg") = 0.);

  py::class_<DisortWrapper>(m, "disort", 
      R"(Wrapper class for the C-DISORT program)")
      .def_readwrite("btemp", &DisortWrapper::btemp, R"(
          Bottom boundary temperature (K), defaults to 0
          )")
      .def_readwrite("ttemp", &DisortWrapper::ttemp, R"(
          Top boundary temperature (K), defaults to 0
          )")
      .def_readwrite("fluor", &DisortWrapper::fluor, R"(
          Intensity of bottom-boundary isotropic illumination, defaults to 0
          )")
      .def_readwrite("albedo", &DisortWrapper::albedo, R"(
          Surface albedo, defaults to 0. Needed if lamber = True
          )")
      .def_readwrite("fisot", &DisortWrapper::fisot, R"( 
          Intensity of top-boundary isotropic illumination, defaults to 0
          )")
      .def_readwrite("fbeam", &DisortWrapper::fbeam, R"(
          Intensity of top-boundary incident parallel beam illumination, defaults to 0
          )")
      .def_readwrite("temis", &DisortWrapper::temis, R"(
          Emissivity of top boundary. Needed if lamber = True, defaults to 0
          )")
      .def_readwrite("umu0", &DisortWrapper::umu0, R"( 
          Cosine of incident beam zenith angle (positive), defaults to 1.
          )")
      .def_readwrite("phi0", &DisortWrapper::phi0, R"( 
          Azimuthal angle of incident beam, defaults to 0.
          )")

      .def(py::init())

      .def("set_header", &DisortWrapper::SetHeader, R"(
          Set header for disort output

          Parameters
          ----------
          arg0 : str
              Header string

          Returns
          -------
          None
          )")

      .def("__str__", &DisortWrapper::ToString)

      .def("set_accuracy", &DisortWrapper::SetAccuracy, R"(
          Set accuracy of disort convergence

          Parameters
          ----------
          arg0 : float
              Accuracy of disort convergence.

          Returns
          -------
          None

          Notes
          -----
          Default accuracy is 1e-6
          )")

      .def("set_flags", &DisortWrapper::SetFlags, R"( 
          Set radiation flags for disort

          Parameters
          ----------
          arg0 : dict
            Dictionary of radiation flags consisting of a key (flag) and a value (True/False).

          Returns
          -------
          DisortWrapper object

          Examples
          --------
          >>> import pydisort
          >>> disort = pydisort.disort()
          >>> dict = {'ibcnd': False, 'usrtau': True, 'usrang': True, 'lamber': True, 'plank': True}
          >>> disort.set_flags(dict).seal()
          >>> rad, flx = disort.run()

          Notes
          -----
          The following flags are supported:

          .. list-table::
             :widths: 25 10 25
             :header-rows: 1

             * - Flag
               - Default value
               - Description
             * - 'ibcnd'
               - False
               - General or Specific boundary condition
             * - 'usrtau'
               - True
               - use user optical depths
             * - 'usrang'
               - True
               - use user azimuthal angles
             * - 'lamber'
               - True
               - turn on lambertian reflection surface
             * - 'plank'
               - True
               - turn on plank source (thermal emission)
             * - 'spher'
               - False
               - turn on spherical correction
             * - 'onlyfl'
               - False
               - only compute radiative fluxes
             * - 'quiet'
               - True
               - turn on disort internal printout
             * - 'intensity_correction'
               - True
               - turn on intensity correction
             * - 'old_intensity_correction'
               - True
               - turn on old intensity correction
             * - 'general_source'
               - False
               - turn on general source
             * - 'output_uum'
               - False
               - output azimuthal components of the intensity
             * - 'print-input'
               - False
               - print input parameters
             * - 'print-fluxes'
               - False
               - print fluxes
             * - 'print-intensity'
               - False
               - print intensity
             * - 'print-transmissivity'
               - False
               - print transmissivity
             * - 'print-phase-function'
               - False
               - print phase function

          A General boundary condition is invoked when 'ibcnd' is set to False.
          This allows:

          - beam illumination from the top (set fbeam)
          - isotropic illumination from the top (set fisot)
          - thermal emission from the top (set ttemp and temis)
          - interal thermal emission (use set_temperature_on_level)
          - reflection at the bottom (set lamber, albedo)
          - thermal emission from the bottom (set btemp)
          
          A Special boundary condition is invoked when 'ibcnd' is set to True. 
          Special boundary condition only returns albedo and transmissivity of 
          the entire medium.

          - current version of pydisort has limited support for this option.
          - consult the documentation of DISORT for more details on this option.
          )")

      .def("set_atmosphere_dimension", &DisortWrapper::SetAtmosphereDimension, R"(
          Set atmosphere dimension for disort

          Parameters
          ----------
          nlyr : int
              Number of layers
          nmom : int, optional, defaults to 4
              Number of phase function moments
          nstr : int, optional, defaults to 4
              Number of streams (computational polar angles)

          Returns
          -------
          DisortWrapper object

          Example
          -------
          >>> import pydisort
          >>> disort = pydisort.disort()
          >>> disort.set_atmosphere_dimension(1, 4, 4, 4)
          >>> print(disort)

          Notes
          -----
          It is common to set nstr = nmom.
          Memory was not allocated until ``disort.seal()`` is called.

          See Also
          --------
          seal : allocate memory for disort
          )",
           py::arg("nlyr"), py::arg("nmom") = 4, py::arg("nstr") = 4)

      .def("set_intensity_dimension", &DisortWrapper::SetIntensityDimension, R"(
          Set intensity dimension for disort

          Parameters
          ----------
          nuphi : int, optional, defaults to 1
              Number of user azimuthal angles
          nutau: int, optional, defaults to 1
              Number of user optical depths
          numu : int, optional, defaults to 1
              Number of user polar angles

          Returns
          -------
          DisortWrapper object

          Example
          -------
          >>> import pydisort
          >>> disort = pydisort.disort()
          >>> disort.set_intensity_dimension(1, 4, 4)

          Notes
          -----
          This function should be called when intensity calculations are requested.
          Memory was not allocated until ``disort.seal()`` is called.

          See Also
          --------
          set_atmosphere_dimension : set atmosphere dimension
          seal : allocate memory for disort
          )",
           py::arg("nuphi") = 1, py::arg("nutau") = 1, py::arg("numu") = 1)

      .def("seal", &DisortWrapper::Seal, R"(
          Seal disort object after setting the dimensions

          Returns
          -------
          None

          Example
          -------
          >>> import pydisort
          >>> disort = pydisort.disort()
          >>> disort.seal()
          >>> print(disort)

          Notes
          -----
          Internal memory is allocated after this function is called.
          )")

      .def("unseal", &DisortWrapper::Unseal, R"(
          Unseal disort object to allow resetting the dimensions

          Returns
          -------
          None

          Example
          -------
          >>> import pydisort
          >>> disort = pydisort.disort()
          >>> disort.seal()
          >>> disort.unseal()
          >>> print(disort)

          Notes
          -----
          Internal memory is deallocated after unseal() is called.
          )")

      .def("set_user_optical_depth", &DisortWrapper::SetUserOpticalDepth, R"(
          Set user optical depth for disort

          Parameters
          ----------
          tau : array_like
              User optical depth

          Returns
          -------
          None
          )")

      .def("set_user_cosine_polar_angle", &DisortWrapper::SetUserCosinePolarAngle, R"(
          Set user cosine polar angle for disort

          Parameters
          ----------
          arg0 : array_like
              User cosine polar angle

          Returns
          -------
          None

          Warning
          -------
          Must be called after ``disort.seal()``
          )")

      .def("set_user_azimuthal_angle", &DisortWrapper::SetUserAzimuthalAngle, R"(
          Set user azimuthal angle for disort

          Parameters
          ----------
          arg0 : array of floats
              User azimuthal angle in radians

          Returns
          -------
          None

          Warning
          -------
          Must be called after ``disort.seal()``
          )")

      .def("set_wavenumber_invcm", &DisortWrapper::SetWavenumber_invcm, R"(
          Set a wavenumber for disort

          Parameters
          ----------
          arg0 : float
              Wavenumber

          Returns
          -------
          None

          Notes
          -----
          This function sets both the minimum and maximum wavenumber to the same value.
          Internal thermal emission is calculated as the planck function at the specified 
          wavenumber.

          See Also
          --------
          set_wavenumber_range_invcm : set a wavenumber range for internal thermal emission
          )")

      .def("set_wavenumber_range_invcm", &DisortWrapper::SetWavenumberRange_invcm,
          R"(
          Set a wavenumber range for disort

          Parameters
          ----------
          wmin : float
              Minimum wavenumber
          wmax : float
              Maximum wavenumber

          Returns
          -------
          None

          Notes
          -----
          This function sets the range of wavenumbers for internal thermal emission.
          Internal thermal emission is calculated as the integral of the planck function
          over the specified wavenumber range.

          See Also
          --------
          set_wavenumber_invcm : set a wavenumber for internal thermal emission
          )",
           py::arg("wmin"), py::arg("wmax"))

      .def("set_optical_thickness", &DisortWrapper::SetOpticalThickness,
          R"(
          Set layer optical thickness

          Parameters
          ----------
          arg0 : array of floats
              Optical thickness

          Returns
          -------
          None

          Example
          -------
          >>> import pydisort
          >>> disort = pydisort.disort()
          >>> disort.set_atmosphere_dimension(4)
          >>> disort.seal()
          >>> disort.set_optical_thickness([1.0, 2.0, 3.0, 4.0])

          Warning
          -------
          This function must be called after ``disort.set_atmosphere_dimension()`` and ``disort.seal()``.

          Notes
          -----
          If the size of the input array is inconsistent with the atmosphere dimension,
          the smaller size is used to set the internal array values.

          See Also
          --------
          set_atmosphere_dimension : set atmosphere dimension
          seal : allocate memory for disort
          )")

      .def("set_single_scattering_albedo", &DisortWrapper::SetSingleScatteringAlbedo,
          R"(
          Set layer single scattering albedo

          Parameters
          ----------
          arg0 : array of floats
              Single scattering albedo

          Returns
          -------
          None

          Example
          -------
          >>> import pydisort
          >>> disort = pydisort.disort()
          >>> disort.set_atmosphere_dimension(4)
          >>> disort.seal()
          >>> disort.set_single_scattering_albedo([0.1, 0.2, 0.3, 0.4])

          Warning
          -------
          This function must be called after ``disort.set_atmosphere_dimension()`` and ``disort.seal()``.

          Notes
          -----
          If the size of the input array is inconsistent with the atmosphere dimension,
          the smaller size is used to set the internal array values.

          See Also
          --------
          set_optical_thickness : set layer optical thickness
          )")

      .def("set_temperature_on_level", &DisortWrapper::SetTemperatureOnLevel,
          R"(
          Set temperature (K) on levels (cell interfaces)

          Parameters
          ----------
          arg0 : array of floats
              Temperature on levels

          Returns
          -------
          None

          Example
          -------
          >>> import pydisort
          >>> disort = pydisort.disort()
          >>> disort.set_atmosphere_dimension(4)
          >>> disort.seal()
          >>> disort.set_temperature_on_level([100.0, 200.0, 300.0, 400.0, 500.0])

          Notes
          -----
          The number of levels is one more than the number of layers.

          Warning
          -------
          This function must be called after ``disort.set_atmosphere_dimension()`` and ``disort.seal()``.

          Notes
          -----
          If the size of the input array is inconsistent with the atmosphere dimension,
          the smaller size is used to set the internal array values.

          See Also
          --------
          set_single_scattering_albedo : set layer single scattering albedo
          )")

      .def("set_phase_moments", 
          [](DisortWrapper &disort, py::array_t<double> &pmom) {
            py::buffer_info info = pmom.request();
            if (info.format != py::format_descriptor<double>::format()) {
              throw std::runtime_error("Incompatible buffer format!");
            } else {
              if (info.ndim == 1) {
                disort.SetPhaseMoments(static_cast<double *>(info.ptr), 1,
                                       info.shape[0]);
              } else if (info.ndim == 2) {
                disort.SetPhaseMoments(static_cast<double *>(info.ptr),
                                       info.shape[0], info.shape[1]);
              } else {
                throw std::runtime_error("Incompatible buffer format!");
              }
            }
          }, R"(
          Set layer phase moments

          Parameters
          ----------
          arg0 : numpy.ndarray, shape (nlyr, nmom + 1) or (nmom + 1,)
              Phase function moments

          Returns
          -------
          None

          Example
          -------
          >>> import pydisort
          >>> ds = pydisort.disort()
          >>> ds.set_atmosphere_dimension(nlyr = 1, nstr = 16, nmom = 16)
          >>> ds.seal()
          >>> _, _, nmom = ds.dimensions()
          >>> pmom = get_phase_function(nmom, "isotropic")
          >>> ds.set_phase_moments(pmom)

          Warning
          -------
          This function must be called after ``disort.set_atmosphere_dimension()`` and ``disort.seal()``.

          Notes
          -----
          If the size of the input array is inconsistent with the atmosphere dimension,
          the smaller size is used to set the internal array values.

          See Also
          --------
          dimensions : get DISORT internal dimensions
          get_phase_function : get phase function moments
          set_single_scattering_albedo : set layer single scattering albedo
          )")

      .def("run",
          [](DisortWrapper &disort) {
            disort.Run();

            py::array_t<double> rad_dummy;

            py::capsule flx_capsule(&disort.Result()->rad[0].rfldir, [](void *) {});
            py::array_t<double> flx({disort.nLayers() + 1, 8}, 
                &disort.Result()->rad[0].rfldir, flx_capsule);

            if (disort.IsFluxOnly()) {
              return std::make_tuple(rad_dummy, flx);
            }

            py::capsule rad_capsule(disort.Result()->uu, [](void *) {});
            py::array_t<double> rad({disort.nUserAzimuthalAngles(), 
                                     disort.nUserOpticalDepths(),
                                     disort.nUserPolarAngles()}, 
                                     disort.Result()->uu,
                                     rad_capsule);
            return std::make_tuple(rad, flx);
          }, R"(
          Run DISORT radiative transfer solver

          Parameters
          ----------
          None

          Returns
          -------
          rad : numpy.ndarray, shape (nuphi, nutau, numu)
              radiation intensity
          flx : numpy.ndarray, shape (nlyr + 1, 8)
              radiation fluxes

          Notes
          -----
          The returned array references the internal memory of the disort object.
          The first dimension is the number of "levels", which is one more than
          the number of layers. The second dimension is the number of flux fields.
          There are 8 flux fields in total, named as follows:
          
          .. list-table::
            :widths: 10 20 30
            :header-rows: 1

            * - Index
              - Name
              - Description
            * - 0
              - RFLDIR
              - Direct beam flux
            * - 1
              - FLDN
              - Diffuse downward flux
            * - 2
              - FLUP
              - Diffuse upward flux
            * - 3
              - DFDT
              - Flux divergence, d (net flux) / d (optical depth)
            * - 4
              - UAVG
              - mean intensity including direct beam
            * - 5
              - UAVGDN
              - mean diffuse downward intensity
            * - 6
              - UAVGUP
              - mean diffuse upward intensity
            * - 7
              - UAVGSO
              - mean direct beam

          All quantities are calculated without delta-M scaling.

          Examples
          --------
          >>> from pydisort import disort
          >>> ds = pydisort.disort()
          >>> ds.seal()
          >>> rad, flx = ds.run()

          In the example above, ``rad`` and ``flx`` references the internal memory of the
          DISORT solver without copying the data. If you want to extract a subset of
          the flux fields, you can do the following:

          >>> from pydisort import disort, RFLDIR, FLDN, FLUP
          >>> ds = pydisort.disort()
          >>> ds.seal()
          >>> _, flx = ds.run()
          >>> flx = flx[:, [RFLDIR, FLDN, FLUP]]

          However, using a subset of the flux fields will create a copy of the
          underlying memory. In the example above, ``flx`` is the variable that extracts 
          the direct beam, diffuse downward, and diffuse upward fluxes.
          )")

      .def("dimensions",
          [](DisortWrapper &disort) {
            std::tuple<int, int, int> shape = {disort.nLayers(), 
                                               disort.nStreams(),
                                               disort.nMoments()};
            return shape;
          }, R"(
          Get the dimensions of the DISORT solver

          Parameters
          ----------
          None

          Returns
          -------
          nlyr : int
              Number of layers
          nstr : int
              Number of streams
          nmom : int
              Number of phase function moments in addition to the zeroth moment
          )");
}
