//
//
//  Created by Administrator on 3/9/2021.
// mex_ECOTRANode_08222022
//  3/4/2021    - changing header file of cubic_b_spline.hpp to cardinal_cubic_b_spline.hpp
//  3/9/2021    - filtering negative values out of biomasses, interpolated physiologies, consumption_IN & predation_OUT, QQQ
//  8/14/2022   - use qb terms to calculate biomasses for physical flux, use pb terms in the dy equation
//              - b) transfer benthic NH4 to sub-surface boxes
//              - c) filter FunctionalResponse. Set negative values to 0
//  9/18/2022   - adding ability for forced input of any group (mmoles N/m3/d), use here to force input NO3 driver


// include h-files******************************************************************
#include <boost/numeric/odeint.hpp> // NOTE: matlab mex compiler doesn't know about boost directory, compile mex with this additional path to headers directory info-->> mex mex_worktest_04182020.cpp -I/usr/local/include/
// #include <boost/math/interpolators/cubic_b_spline.hpp> // NOTE: matlab mex compiler doesn't know about boost directory, compile mex with this additional path to headers directory info-->> mex mex_INTERPOLATIONtest_03312020.cpp -I/usr/local/include/
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp> // QQQ NOTE: matlab mex compiler doesn't know about boost directory, compile mex with this additional path to headers directory info-->> mex mex_INTERPOLATIONtest_03312020.cpp -I/usr/local/include/
#include "mex.hpp"
#include "mexAdapter.hpp"
// *********************************************************************************




// declarations*********************************************************************
// define the type of container used to hold the ODE state as a vector variable type using the vector template within the Standard Library (std); other state_types than vector are possible, boost::array<double, N_dimms>
typedef std::vector<double> state_type;

// define odeint stepper
boost::numeric::odeint::runge_kutta_dopri5<state_type> stepper;

// Get array factory to return solution to MATLAB
matlab::data::ArrayFactory factory; // Use ArrayFactory to create matlab::data::Array objects.
// *********************************************************************************





// SplineInterpolate class**********************************************************
class SplineInterpolate {
    // interpolate time-dependent terms @ t
    // NOTE: currently only works with 3-D matrices, where time is given along the row axis & num_rows must = 1, num_t is > 1
    // NOTE: usage to instantiate a physiological parameter time-series to interpolate: SplineInterpolate METABOLISM(timeseries, num_t, num_rows, num_clms, num_layers, step);
    // includes methods:
    //          none
    // FFF can make work with 4D matrices in future
    
private:
    int c_num_t, c_num_rows, c_num_clms, c_num_layers;
    double step;
    std::vector<double>& timeseries;
    double* current_timeseries; // create pointer variable "current_timeseries" of type double* that points to dynamically allocated array
    
public:
    // boost::math::cubic_b_spline<double>* spline; // create pointer variable "spline" of type boost::math::cubic_b_spline<double>* that points to dynamically allocated array
    boost::math::interpolators::cardinal_cubic_b_spline<double>* spline; // QQQ create pointer variable "spline" of type boost::math::cubic_b_spline<double>* that points to dynamically allocated array

    
    // constructor -----------------------------------------------------------------
    SplineInterpolate(std::vector<double>& timeseries_in, int num_t_in, int num_rows_in, int num_clms_in, int num_layers_in, double dt_in) :
        timeseries{timeseries_in},
        step{dt_in},
        c_num_t{num_t_in},
        c_num_rows{num_rows_in},
        c_num_clms{num_clms_in},
        c_num_layers{num_layers_in} {
        
        current_timeseries  = new double[c_num_t]; // allocate num_t type double elements on heap memory
        // spline              = new boost::math::cubic_b_spline<double>[c_num_rows * c_num_clms * c_num_layers]; // allocate (num_rows * num_clms * num_layers) type boost::math::cubic_b_spline<double> elements on heap
        spline              = new boost::math::interpolators::cardinal_cubic_b_spline<double>[c_num_rows * c_num_clms * c_num_layers]; // QQQ allocate (num_rows * num_clms * num_layers) type boost::math::cubic_b_spline<double> elements on heap

        
        for (int layer_index=0; layer_index<(c_num_layers); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_clms); ++clm_index) {
                for (int row_index=0; row_index<(c_num_rows); ++row_index) {
                    
                    for (int t_index=0; t_index<(c_num_t); ++t_index) {
                        current_timeseries[t_index] = timeseries[(layer_index * c_num_t * c_num_clms) + (t_index * c_num_clms) + clm_index];
                    } // end for t_index
                    
                    // spline[(layer_index * c_num_rows * c_num_clms) + (row_index * c_num_clms) + clm_index] = boost::math::cubic_b_spline<double>(current_timeseries, c_num_t, 1 /* start time */, step); // NOTE: start time t = 1 (because t_grid(1) = 1 from MATLAB)
                    spline[(layer_index * c_num_rows * c_num_clms) + (row_index * c_num_clms) + clm_index] = boost::math::interpolators::cardinal_cubic_b_spline<double>(current_timeseries, c_num_t, 1 /* start time */, step); // QQQ NOTE: start time t = 1 (because t_grid(1) = 1 from MATLAB)

                    
                } // end for row_index
            } // end for clm_index
        } // end for layer_index
            
    std::cout << "Constructed class SplineInterpolate" << std::endl;
    } // end constructor -----------------------------------------------------------
    
    
    // destructor ------------------------------------------------------------------
    ~SplineInterpolate() {
        delete[] spline;
        delete[] current_timeseries;
        std::cout << "Destroyed class SplineInterpolate" << std::endl;
    } // NOTE: destructor is automatically called for each instantiated object at end of main()
    
}; // end class SplineInterpolate
// *********************************************************************************





// PhysicalFlux class***************************************************************
class PhysicalFlux {
    // includes methods:
    //          1) m_PhysicalFlux_intraODE
    
private:
    std::size_t array_dims;
    unsigned int index_a, index_b, index_c, index_d, index_looky;
    double* cp_compact_biomass_source;
    double* cp_biomass_boundary_import;
    // double* cp_biomass_boundary_export; // option for domain export values
    double* cp_flux_biomass;
    double* cp_flux_domain_import;
    // double* cp_flux_domain_export; // option for domain export values

public:
    // shared variables, constant over time
    static unsigned int cs_num_grps, cs_num_boxes, cs_num_drivers;
    static std::vector<unsigned int> cs_looky_driver; // row address(es) of driver group(s) (e.g., NO3); (vector: num_drivers X 1); QQQ is this declaration necessary if I use PhysicaFlux::cs_looky_driver ?
    static std::vector<unsigned int> cs_looky_externalForcing; // row address(es) of external forcing group(s) (e.g., NO3, juv salmon); (vector: num_externalForcing_grps X 1); QQQ is this declaration necessary if I use PhysicaFlux::cs_looky_driver ?

    // shared variables, time-dependent
    static std::vector<double> cs_biomass_t; // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
    static std::vector<double> cs_BoxVolume_t; // (m3); (vector: 1 X num_boxes)
    
    // flux-specific variables, constant over time
    unsigned int c_num_fluxes, c_num_fluxes_BoundaryImport; // (NOTE: const does not work for these two variables, constructor needs explicitly defined values when const) QQQ make private?
    // unsigned int c_num_fluxes_BoundaryExport; // option for domain export values
    std::vector<unsigned int> c_looky_flux; // (2D matrix: num_fluxes X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links QQQ make private?
    std::vector<unsigned int> c_looky_boundary_import; // (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
    // std::vector<unsigned int> c_looky_boundary_export; // option for domain export values; (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
    std::vector<double> c_RetentionScaler; // (2D matrix: num_grps X num_boxes DESTINY); NOTE: c_RetentionScaler is kept as all 0 for sinking and migration fluxes
    
    // flux-specific variables, time-dependent
    std::vector<double> c_FLUX_compact_t; // (m3/d); (2D matrix: num_grps X num_fluxes); NOTE: num_fluxes = number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
    std::vector<double> c_flux_import_t; // pre-initialized as 0 each time-step; (mmoles N/d) & (mmoles N/m3/d); (2D matrix: num_grps X (num_boxes+1))
    std::vector<double> c_flux_export_t; // pre-initialized as 0 each time-step; (mmoles N/d) & (mmoles N/m3/d); (2D matrix: num_grps X (num_boxes+1))
    std::vector<double> c_flux_domain_import_t; // pre-initialized as 0 each time-step; (2D matrix: num_grps X (num_boxes+1))
    // std::vector<double> c_flux_domain_export_t; // option for domain export values; pre-initialized as 0 each time-step; (2D matrix: num_grps X (num_boxes+1))
    std::vector<double> c_flux_domain_import_driver_t; // pre-initialized as 0 each time-step; (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_flux_net_t; // pre-initialized as 0 each time-step; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_biomass_boundary_t; // boundary biomass conditions @ t; (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1)); NOTE: c_biomass_boundary_t is set to all 0 for sinking fluxes and is set independenntly for migration fluxes

    
    
    // constructor -----------------------------------------------------------------
    PhysicalFlux(unsigned int num_fluxes_in, unsigned int num_fluxes_BoundaryImport_in /* , unsigned int num_fluxes_BoundaryExport_in */) :
        c_num_fluxes{num_fluxes_in},
        c_num_fluxes_BoundaryImport{num_fluxes_BoundaryImport_in} /*,
        c_num_fluxes_BoundaryExport{num_fluxes_BoundaryExport_in} // NOTE option for domain export values */ {

        cp_compact_biomass_source  = new double[cs_num_grps * c_num_fluxes]{0}; // biomasses of each source box exporting water; (mmoles N/m3); (2D matrix: num_grps X num_fluxes)
        cp_biomass_boundary_import = new double[cs_num_grps * c_num_fluxes_BoundaryImport]{0}; // biomasses immediately outside the model domain of each destiny box that imports water from outside the model domain; (mmoles N/m3); (2D matrix: num_grps X num_fluxes_BoundaryImport); sinking units (mmoles N/m2/d)
        // cp_biomass_boundary_export = new double[cs_num_grps * c_num_fluxes_BoundaryExport]{0}; // option for domain export values; biomasses of each source box exporting water out of domain; (mmoles N/m3); (2D matrix: num_grps X num_export_fluxes); sinking units (mmoles N/m2/d)
        cp_flux_biomass            = new double[cs_num_grps * c_num_fluxes]{0}; // flux biomasses of each source box exporting water; (mmoles N/d); (2D matrix: num_grps X num_fluxes); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d)
        cp_flux_domain_import      = new double[cs_num_grps * c_num_fluxes_BoundaryImport]{0}; // flux biomasses into each destiny box that imports water from outside the model domain; (mmoles N/d); (2D matrix: num_grps X num_fluxes_BoundaryImport); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d)
        // cp_flux_domain_export      = new double[cs_num_grps * c_num_fluxes_BoundaryExport]{0}; // flux biomasses out of each source box that exports water out of the model domain; (mmoles N/d); (2D matrix: cs_num_grps X c_num_fluxes_BoundaryExport); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d); option for domain export values
        
        // initialize as 0 at correct (and unchanging) matrix sizes
        c_looky_flux.resize(c_num_fluxes * 3); // (2D matrix: num_fluxes X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links
        c_looky_boundary_import.resize(c_num_fluxes_BoundaryImport * 3); // (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        // c_looky_boundary_export.resize(c_num_fluxes_BoundaryExport * 3); // option for domain export values; (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        c_RetentionScaler.resize(cs_num_grps * cs_num_boxes); // initialize vector size to all 0; (2D matrix: cs_num_grps X cs_num_boxes DESTINY)
        c_biomass_boundary_t.resize(cs_num_grps * 1 * (cs_num_boxes+1)); // initialize vector size to all 0; boundary biomass conditions @ t; (mmoles N/m3); (3D matrix: cs_num_grps X 1 X (cs_num_boxes+1))
        c_FLUX_compact_t.resize(cs_num_grps * c_num_fluxes); // (m3/d); (horizontal vector: cs_num_grps X c_num_fluxes)
        c_flux_import_t.resize(cs_num_grps * (cs_num_boxes+1)); // pre-initialized as 0 each time-step; (mmoles N/d) & (mmoles N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        c_flux_export_t.resize(cs_num_grps * (cs_num_boxes+1)); // pre-initialized as 0 each time-step; (mmoles N/d) & (mmoles N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        c_flux_domain_import_t.resize(cs_num_grps * (cs_num_boxes+1)); // pre-initialized as 0 each time-step; (2D matrix: num_grps X (num_boxes+1))
        // c_flux_domain_export_t.resize(cs_num_grps * (cs_num_boxes+1)); // option for domain export values; pre-initialized as 0 each time-step; (2D matrix: num_grps X (num_boxes+1))
        c_flux_domain_import_driver_t.resize(1 * cs_num_grps * cs_num_boxes); // pre-initialized as 0 each time-step; (3D matrix: 1 X num_grps X num_boxes)
        c_flux_net_t.resize(1 * cs_num_grps * cs_num_boxes); // pre-initialized as 0 each time-step; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
        
        std::cout << "Constructed class PhysicalFlux" << std::endl;

    } // end constructor -----------------------------------------------------------

    
    
    // destructor ------------------------------------------------------------------
    ~PhysicalFlux() {
        delete[] cp_compact_biomass_source;
        delete[] cp_biomass_boundary_import;
        // delete[] cp_biomass_boundary_export; // option for domain export values
        delete[] cp_flux_biomass;
        delete[] cp_flux_domain_import;
        // delete[] cp_flux_domain_export; // option for domain export values
        
        std::cout << "Destroyed class PhysicalFlux" << std::endl;
    } // NOTE: destructor is automatically called for each instantiated object at end of main()

    
    
    // method 1: m_PhysicalFlux_intraODE ---------------------------------------------
    void m_PhysicalFlux_intraODE() {
//        std::cout << "Running PhysicalFlux.m_PhysicalFlux_intraODE" << std::endl;
        
        for (int clm_index=0; clm_index<(c_num_fluxes); ++clm_index) {
            for (int row_index=0; row_index<(cs_num_grps); ++row_index) {
                
                // STEP 1: prepare box biomass terms--------------------------------
                //         intra-domain source & export box biomasses
                //         3D matrix indexing is: [(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
                index_a     = (0 * cs_num_grps * c_num_fluxes) + (row_index * c_num_fluxes) + clm_index; // index into cp_compact_biomass_source (2D matrix: num_grps X num_fluxes)
                index_b     = (0 * c_num_fluxes * 3) + (clm_index * 3) + 1; // index into c_looky_flux(:, 2) (2D matrix: num_fluxes X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)])
                index_looky = c_looky_flux[index_b];
                index_c     = (index_looky * cs_num_grps * 1) + (row_index * 1) + 0; // index into cs_biomass_t (3D matrix: num_grps X 1 X (num_boxes+1))
                
                cp_compact_biomass_source[index_a] = cs_biomass_t[index_c]; // biomasses of each source box exporting water; (mmoles N/m3); (2D matrix: num_grps X num_fluxes); sinking units (mmoles N/m2/d)
                // -----------------------------------------------------------------

    
                // STEP 2: express fluxes in terms of biomass exchange--------------
                index_a = (0 * cs_num_grps * c_num_fluxes) + (row_index * c_num_fluxes) + clm_index; // index into cp_compact_biomass_source, cp_flux_biomass, & c_FLUX_compact_t (2D matrix: num_grps X num_fluxes)

                cp_flux_biomass[index_a] = c_FLUX_compact_t[index_a] * cp_compact_biomass_source[index_a]; // flux biomasses of each source box exporting water; (mmoles N/d); (2D matrix: num_grps X num_fluxes); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d)
                // -----------------------------------------------------------------

    
                // STEP 3: pool fluxes in SOURCE & DESTINY boxes--------------------
                //         NOTE: this code performs the role of MATLAB accumarray
                index_a     = (0 * cs_num_grps * c_num_fluxes) + (row_index * c_num_fluxes) + clm_index; // index into cp_flux_biomass (2D matrix: num_grps X num_fluxes)
                
                // pool import fluxes
                index_b     = (0 * c_num_fluxes * 3) + (clm_index * 3) + 0; // index into looky_flux DESTINY box (clm = 0) (2D matrix: num_fluxes X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links
                index_looky = c_looky_flux[index_b];
                index_c     = (0 * cs_num_grps * (cs_num_boxes+1)) + (row_index * (cs_num_boxes+1)) + index_looky; // index into flux_import_t (2D matrix: num_grps X (num_boxes+1))
                
                c_flux_import_t[index_c] = c_flux_import_t[index_c] + cp_flux_biomass[index_a]; // pool import to destiny boxes; (mmoles N/d); (3D matrix: num_grps X (num_boxes+1)); NOTE: does NOT include domain import (see cp_flux_domain_import for domain input rates)
                
                // pool export fluxes
                index_b     = (0 * c_num_fluxes * 3) + (clm_index * 3) + 1; // index into looky_flux SOURCE box (clm = 1) (2D matrix: num_fluxes X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links
                index_looky = c_looky_flux[index_b];
                index_c     = (0 * cs_num_grps * (cs_num_boxes+1)) + (row_index * (cs_num_boxes+1)) + index_looky; // index into flux_export_t (2D matrix: num_grps X (num_boxes+1))
                
                c_flux_export_t[index_c] = c_flux_export_t[index_c] + cp_flux_biomass[index_a]; // pool export from source boxes; (mmoles N/d); (2D matrix: num_grps X (num_boxes+1)); NOTE: does include domain export
                // -----------------------------------------------------------------
                
            } // end for row_index
        } // end for clm_index

    
        for (int clm_index=0; clm_index<(c_num_fluxes_BoundaryImport); ++clm_index) {
            for (int row_index=0; row_index<(cs_num_grps); ++row_index) {

                // STEP 4: boundary conditions for import boxes --------------------
                index_a     = (0 * cs_num_grps * c_num_fluxes_BoundaryImport) + (row_index * c_num_fluxes_BoundaryImport) + clm_index; // index into cp_biomass_boundary_import (2D matrix: num_grps X num_fluxes_BoundaryImport)
                index_b     = (0 * c_num_fluxes_BoundaryImport * 3) + (clm_index * 3) + 0; // index into looky_boundary_import (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
                index_looky = c_looky_boundary_import[index_b];
                index_c     = (index_looky * cs_num_grps * 1) + (row_index * 1) + 0; // index into c_biomass_boundary_t (3D matrix: num_grps X 1 X (num_boxes+1))

                cp_biomass_boundary_import[index_a] = c_biomass_boundary_t[index_c]; // biomasses immediately outside the model domain of each destiny box that imports water from outside the model domain; (mmoles N/m3); (2D matrix: num_grps X num_fluxes_BoundaryImport); sinking units (mmoles N/m2/d)
                // -----------------------------------------------------------------


                // STEP 5: express fluxes in terms of biomass exchange--------------
                index_a     = (0 * c_num_fluxes_BoundaryImport * 3) + (clm_index * 3) + 2; // index into looky_boundary_import (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: always drawing from the third column (addresses of import flux clms in compact_flux_t)
                index_looky = c_looky_boundary_import[index_a];
                index_b     = (0 * cs_num_grps * c_num_fluxes) + (row_index * c_num_fluxes) + index_looky; // index into c_FLUX_compact_t (2D matrix: num_grps X num_fluxes)
                index_c     = (0 * cs_num_grps * c_num_fluxes_BoundaryImport) + (row_index * c_num_fluxes_BoundaryImport) + clm_index; // index into cp_biomass_boundary_import & cp_flux_domain_import (2D matrix: num_grps X num_fluxes_BoundaryImport)

                cp_flux_domain_import[index_c] = c_FLUX_compact_t[index_b] * cp_biomass_boundary_import[index_c]; // flux biomasses into each destiny box that imports water from outside the model domain; (mmoles N/d); (2D matrix: num_grps X num_fluxes_BoundaryImport); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d)
                // -----------------------------------------------------------------


                // STEP 6: prepare cross-domain boundary fluxes---------------------
                index_a     = (0 * c_num_fluxes_BoundaryImport * 3) + (clm_index * 3) + 0; // index into looky_boundary_import (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
                index_looky = c_looky_boundary_import[index_a];
                index_b     = (0 * cs_num_grps * c_num_fluxes_BoundaryImport) + (row_index * c_num_fluxes_BoundaryImport) + clm_index; // index into cp_flux_domain_import (2D matrix: num_grps X num_fluxes_BoundaryImport)
                index_c     = (0 * cs_num_grps * (cs_num_boxes+1)) + (row_index * (cs_num_boxes+1)) + index_looky; // index into flux_domain_import_t (2D matrix: num_grps X num_boxes+1)

                c_flux_domain_import_t[index_c] = cp_flux_domain_import[index_b]; // (mmoles N/d); (2D matrix: num_grps X (num_boxes+1))
                // -----------------------------------------------------------------

            } // end for row_index
        } // end for clm_index

        
        // option for domain export values
        //        for (int clm_index=0; clm_index<(c_num_fluxes_BoundaryExport); ++clm_index) {
        //            for (int row_index=0; row_index<(cs_num_grps); ++row_index) {
        //                // 3D matrix indexing is: [(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        //
        //                // STEP 7: boundary conditions for export boxes --------------------
        //                index_a     = (0 * cs_num_grps * c_num_fluxes_BoundaryExport) + (row_index * c_num_fluxes_BoundaryExport) + clm_index; // index into cp_biomass_boundary_export (2D matrix: cs_num_grps X c_num_fluxes_BoundaryExport)
        //                index_b     = (0 * c_num_fluxes_BoundaryExport * 3) + (clm_index * 3) + 1; // index into c_looky_boundary_export (2D matrix: c_num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)]); NOTE: always drawing from the second column (SOURCE)
        //                index_looky = c_looky_boundary_export[index_b];
        //                index_c     = (index_looky * cs_num_grps * 1) + (row_index * 1) + 0; // index into cs_biomass_t (3D matrix: num_grps X 1 X (num_boxes+1))
        //
        //                cp_biomass_boundary_export[index_a] = cs_biomass_t[index_c]; // biomasses of each source box exporting water out of domain; (mmoles N/m3); (2D matrix: num_grps X c_num_fluxes_BoundaryExport); sinking units (mmoles N/m2/d)
        //                // -----------------------------------------------------------------
        //
        //
        //                // STEP 8: express fluxes in terms of biomass exchange--------------
        //                index_a     = (0 * c_num_fluxes_BoundaryExport * 3) + (clm_index * 3) + 2; // index into c_looky_boundary_export (2D matrix: c_num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)]); NOTE: always drawing from the third column (addresses of import flux clms in compact_flux_t)
        //                index_looky = c_looky_boundary_export[index_a];
        //                index_b     = (0 * cs_num_grps * c_num_fluxes) + (row_index * c_num_fluxes) + index_looky; // index into c_FLUX_compact_t (2D matrix: num_grps X num_fluxes)
        //                index_c     = (0 * cs_num_grps * c_num_fluxes_BoundaryExport) + (row_index * c_num_fluxes_BoundaryExport) + clm_index; // index into cp_biomass_boundary_export & cp_flux_domain_export (2D matrix: cs_num_grps X c_num_fluxes_BoundaryExport)
        //
        //                cp_flux_domain_export[index_c] = c_FLUX_compact_t[index_b] * cp_biomass_boundary_export[index_c]; // flux biomasses out of each source box that exports water out of the model domain; (mmoles N/d); (2D matrix: cs_num_grps X c_num_fluxes_BoundaryExport); sinking units (m2)*(mmoles N/m2/d) = (mmoles N/d)
        //                // -----------------------------------------------------------------
        //
        //
        //                // STEP 9: prepare cross-domain boundary fluxes---------------------
        //                // 3D matrix indexing is: [(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        //                index_a     = (0 * c_num_fluxes_BoundaryExport * 3) + (clm_index * 3) + 1; // index into c_looky_boundary_export (2D matrix: c_num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)]); NOTE: always drawing from the second column (SOURCE)
        //                index_looky = c_looky_boundary_export[index_a];
        //                index_b     = (0 * cs_num_grps * c_num_fluxes_BoundaryExport) + (row_index * c_num_fluxes_BoundaryExport) + clm_index; // index into cp_flux_domain_export (2D matrix: cs_num_grps X c_num_fluxes_BoundaryExport)
        //                index_c     = (0 * cs_num_grps * (cs_num_boxes+1)) + (row_index * (cs_num_boxes+1)) + index_looky; // index into c_flux_domain_export_t (2D matrix: num_grps X (num_boxes+1))
        //
        //                c_flux_domain_export_t[index_c] = cp_flux_domain_export[index_b]; // (mmoles N/d); (2D matrix: num_grps X (num_boxes+1))
        //                // -----------------------------------------------------------------
        //
        //            } // end for row_index
        //        } // end for clm_index


        // STEP 10: calculate driver group production (flux_domain_import_driver_t)
        //         NOTE: keep production_input groups (e.g., NO3) separate from other
        //                boundary condition groups as driver_input_t is treated
        //                separately in dy calculation
        //         3D matrix indexing is: [(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        for (int row_index=0; row_index<(cs_num_drivers); ++row_index) {

            index_a = (0 * cs_num_drivers * 1) + (row_index * 1) + 0; // index into looky_driver (num_drivers)
            index_looky = cs_looky_driver[index_a]; // QQQ should this be PhysicalFlux::cs_looky_driver[index_a] ?

            for (int clm_index=0; clm_index<(cs_num_boxes); ++clm_index) {

                index_b = (clm_index * 1 * cs_num_grps) + (0 * cs_num_grps) + index_looky; // index into flux_domain_import_driver_t (3D matrix: 1 X num_grps X num_boxes)
                index_c = (0 * cs_num_grps * (cs_num_boxes+1)) + (index_looky * (cs_num_boxes+1)) + clm_index; // index into flux_domain_import_t (2D matrix: num_grps X (num_boxes+1) DESTINY) NOTE: clm_index skips over the box+1 boundary border column

                c_flux_domain_import_driver_t[index_b]  = c_flux_domain_import_t[index_c]; // (mmoles N/d); (3D matrix: 1 X num_grps X num_boxes DESTINY)
                c_flux_domain_import_t[index_c]         = 0; // zero out driver group(s); (mmoles N/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)

            } // end for clm_index (NOTE order swap between row and clm loops)
        } // end for row_index
        // -------------------------------------------------------------------------


        // STEP 11: apply retention scalers, calculate flux relative to destiny box volume
        // 3D matrix indexing is: [(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        for (int clm_index=0; clm_index<(cs_num_boxes); ++clm_index) {
            for (int row_index=0; row_index<(cs_num_grps); ++row_index) {

                index_a = (0 * cs_num_grps * (cs_num_boxes+1)) + (row_index * (cs_num_boxes+1)) + clm_index; // index into flux_import_t, flux_export_t, flux_domain_import_t (2D matrix: num_grps X (num_boxes+1) DESTINY); NOTE: clm_index does not consider the last (+1 boundary) column (in MATLAB, the last column is deleted)
                index_b = (0 * cs_num_grps * cs_num_boxes) + (row_index * cs_num_boxes) + clm_index; // index into RetentionScaler (2D matrix: num_grps X num_boxes DESTINY)
                index_c = (0 * 1 * cs_num_boxes) + (0 * cs_num_boxes) + clm_index; // index into BoxVolume_t (vector: 1 X num_boxes)
                index_d = (clm_index * 1 * cs_num_grps) + (0 * cs_num_grps) + row_index; // index into flux_net_t & flux_domain_import_driver_t (3D matrix: 1 X num_grps X num_boxes)

                // apply retention scalers
                c_flux_import_t[index_a]                  = c_flux_import_t[index_a]               * (1 - c_RetentionScaler[index_b]); // (mmoles N/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)
                c_flux_export_t[index_a]                  = c_flux_export_t[index_a]               * (1 - c_RetentionScaler[index_b]); // (mmoles N/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)
                c_flux_domain_import_t[index_a]           = c_flux_domain_import_t[index_a]        * (1 - c_RetentionScaler[index_b]); // (mmoles N/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)
                // c_flux_domain_export_t[index_a]           = c_flux_domain_export_t[index_a]        * (1 - c_RetentionScaler[index_b]); // option for domain export values; (mmoles N/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)
                c_flux_domain_import_driver_t[index_d]    = c_flux_domain_import_driver_t[index_d] * (1 - c_RetentionScaler[index_b]); // (mmoles N/d); (3D matrix: 1 X num_grps X num_boxes)
                
                // calculate flux relative to destiny box volume
                c_flux_import_t[index_a]                  = c_flux_import_t[index_a]               / cs_BoxVolume_t[index_c]; // (mmoles N/m3/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)
                c_flux_export_t[index_a]                  = c_flux_export_t[index_a]               / cs_BoxVolume_t[index_c]; // (mmoles N/m3/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)
                c_flux_domain_import_t[index_a]           = c_flux_domain_import_t[index_a]        / cs_BoxVolume_t[index_c]; // (mmoles N/m3/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)
                // c_flux_domain_export_t[index_a]           = c_flux_domain_export_t[index_a]        / cs_BoxVolume_t[index_c]; // option for domain export values; (mmoles N/m3/d); (2D matrix: num_grps X (num_boxes+1) DESTINY)
                c_flux_domain_import_driver_t[index_d]    = c_flux_domain_import_driver_t[index_d] / cs_BoxVolume_t[index_c]; // (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes)

                // calculate flux_net_t
                c_flux_net_t[index_d] = (c_flux_import_t[index_a] + c_flux_domain_import_t[index_a]) - c_flux_export_t[index_a]; // net biomass import to each box; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)

            } // end for row_index
        } // end for clm_index
        // -------------------------------------------------------------------------
        
    } // end method 1: m_PhysicalFlux_intraODE--------------------------------------

}; // end class physical flux
// *********************************************************************************





// ECOTRANode class*****************************************************************
class ECOTRANode {
    // includes 6 methods:
    //              1) m_interpolate
    //              2) m_calc_Biomass_t
    //              3) m_calc_PhysicalFlux_t
    //              4) m_calc_EnergyBudget_t
    //              5) m_calc_CONSUMPTION_t
    //              6) (m_calc_dxdt_t) uses overloaded operator()

private:
    std::size_t array_dims;
    unsigned int index_a, index_b, index_c, index_looky, index_fate, index_ConsumptionBudget, index_EnergyBudget, index_looky_plgcNH4, index_looky_bnthNH4;
    double temp_biomass, feces_t, senescence_t, numerator, denominator, temp_sum_row, temp_sum_clm, temp_sum_bnthNH4;
    
public:
    unsigned int c_num_grps, c_num_boxes, c_num_drivers, c_num_externalForcing_grps; // NEW!!!
    unsigned char c_num_nutrients, c_num_eggs, c_num_ANYdetritus, c_num_livingANDfleets, c_num_bnthNH4;
    
    double t;
    
    std::vector<unsigned char> c_looky_nutrients, c_looky_eggs, c_looky_livingANDfleets, c_looky_ANYdetritus, c_looky_plgcNH4, c_looky_bnthNH4; // (vector: 1 X num_QQQ)
    
    std::vector<double> c_FunctionalResponseParams; // (3D matrix: num_grps X num_grps X num_boxes)
    
    std::vector<double> c_production_initial; // (3D matrix: num_grps X 1 X num_boxes); NOTE: used to define initial state conditions and as baseline conditions for functional response calculations
    
    std::vector<double> c_TransferEfficiency; // (3D matrix: 1 X num_grps X num_boxes)
    
    //    state_type dxdt; // (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes); only needed for debugging, we don't need this when ECOTRANode instance is being called by ODEint

    
    std::vector<double> c_EnergyBudget_t; // (3D matrix: num_grps X num_grps X num_boxes)

    std::vector<double> c_fate_feces; // (3D matrix: num_ANYdetritus X num_grps X num_boxes)
    std::vector<double> c_fate_metabolism; // (3D matrix: num_nutrients X num_grps X num_boxes)
    std::vector<double> c_fate_eggs; // (3D matrix: num_eggs X num_grps X num_boxes)
    std::vector<double> c_fate_predation; // (3D matrix: num_livingANDfleets X num_grps X num_boxes)
    std::vector<double> c_fate_senescence; // (3D matrix: num_ANYdetritus X num_grps X num_boxes)


    std::vector<double> c2_biomass_boundary_t;          // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1)); declare 1 block of dynamic memory
    std::vector<double> c_biomass_MigratorBoundary_t;   // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1)); declare 1 block of dynamic memory
    std::vector<double> c_external_driver_t;            // (3D matrix: 1 X num_drivers X (num_boxes+1))
    std::vector<double> c_externalForcing_grps_t;       // NEW!!! (3D matrix: 1 X num_externalForcing_grps X num_boxes)
    std::vector<double> c_externalForcing_t;            // NEW!!! pre-initialized as 0 each time-step; (3D matrix: 1 X num_grps X num_boxes)


    std::vector<double> c_ConsumptionBudget_feces_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_ConsumptionBudget_metabolism_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_ConsumptionBudget_eggs_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_ConsumptionBudget_predation_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_ConsumptionBudget_senescence_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_ConsumptionBudget_ba_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_ConsumptionBudget_em_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_qb_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_pb_t; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_ThorntonLessemScaler_t;         // (3D matrix: 1 X num_grps X num_boxes)

    std::vector<double> consumption_IN; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> predation_OUT; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> dy; // (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> flux_domain_import_driver_total_t; // (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
    std::vector<double> c_FunctionalResponse; // (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
    std::vector<double> Q_cp_t; // (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)

    // pointers to SplineInterpolate class & PhysicalFlux class
    std::vector<SplineInterpolate*> c_ptr_InterpolatePhysiology;
    std::vector<SplineInterpolate*> c_ptr_InterpolatePhysicalFlux;
    std::vector<PhysicalFlux*> c_ptr_PhysicalFlux;

    std::vector<double> c_preinit_flux_import_t; // (2D matrix: num_grps X (num_boxes+1))
    std::vector<double> c_preinit_flux_export_t; // (2D matrix: num_grps X (num_boxes+1))
    std::vector<double> c_preinit_flux_domain_import_t; // (2D matrix: num_grps X (num_boxes+1))
    // std::vector<double> c_preinit_flux_domain_export_t; // option for domain export values; (2D matrix: num_grps X (num_boxes+1))
    std::vector<double> c_preinit_flux_domain_import_driver_t; // (3D matrix: 1 X num_grps X num_boxes DESTINY)
    std::vector<double> c_preinit_flux_net_t; // (3D matrix: 1 X num_grps X num_boxes DESTINY)

    std::vector<unsigned int> c_num_fluxes; // (vector: 11 X 1)
    
    
    // constructor -----------------------------------------------------------------
    ECOTRANode(unsigned int num_grps_in,
               unsigned char num_nutrients_in,
               unsigned char num_bnthNH4_in,
               unsigned char num_eggs_in,
               unsigned char num_ANYdetritus_in,
               unsigned char num_livingANDfleets_in,
               unsigned int num_boxes_in,
               unsigned int num_drivers_in,
               unsigned int num_externalForcing_grps_in,
               std::vector<unsigned int> num_fluxes_in,
               std::vector<unsigned char> looky_nutrients_in,
               std::vector<unsigned char> looky_eggs_in,
               std::vector<unsigned char> looky_livingANDfleets_in,
               std::vector<unsigned char> looky_ANYdetritus_in,
               std::vector<unsigned char> looky_plgcNH4_in, // QQQ what if num_plgcNH4 = 0???; (vector: 1 X num_plgcNH4); FFF replace with looky_NH4
               std::vector<unsigned char> looky_bnthNH4_in, // QQQ what if num_bnthNH4 = 0???; (vector: 1 X num_bnthNH4); FFF replace with looky_NH4
               std::vector<SplineInterpolate*> ptr_InterpolatePhysiology_in,
               std::vector<SplineInterpolate*> ptr_InterpolatePhysicalFlux_in,
               std::vector<PhysicalFlux*> ptr_PhysicalFlux_in,
               std::vector<double> EnergyBudget_in,
               std::vector<double> fate_feces_in,
               std::vector<double> fate_metabolism_in,
               std::vector<double> fate_eggs_in,
               std::vector<double> fate_predation_in,
               std::vector<double> fate_senescence_in,
               std::vector<double> TransferEfficiency_in,
               std::vector<double> FunctionalResponseParams_in,
               std::vector<double> production_initial_in) :
            c_num_grps{num_grps_in},
            c_num_nutrients{num_nutrients_in},
            c_num_bnthNH4{num_bnthNH4_in},
            c_num_eggs{num_eggs_in},
            c_num_ANYdetritus{num_ANYdetritus_in},
            c_num_livingANDfleets{num_livingANDfleets_in},
            c_num_boxes{num_boxes_in},
            c_num_drivers{num_drivers_in},
            c_num_externalForcing_grps{num_externalForcing_grps_in},
            c_num_fluxes{num_fluxes_in},
            c_looky_nutrients{looky_nutrients_in},
            c_looky_eggs{looky_eggs_in},
            c_looky_livingANDfleets{looky_livingANDfleets_in},
            c_looky_ANYdetritus{looky_ANYdetritus_in},
            c_looky_plgcNH4{looky_plgcNH4_in},
            c_looky_bnthNH4{looky_bnthNH4_in},
            c_ptr_InterpolatePhysiology{ptr_InterpolatePhysiology_in},
            c_ptr_InterpolatePhysicalFlux{ptr_InterpolatePhysicalFlux_in},
            c_ptr_PhysicalFlux{ptr_PhysicalFlux_in},
            c_EnergyBudget_t{EnergyBudget_in},
            c_fate_feces{fate_feces_in},
            c_fate_metabolism{fate_metabolism_in},
            c_fate_eggs{fate_eggs_in},
            c_fate_predation{fate_predation_in},
            c_fate_senescence{fate_senescence_in},
            c_TransferEfficiency{TransferEfficiency_in},
            c_FunctionalResponseParams{FunctionalResponseParams_in},
            c_production_initial{production_initial_in} {
            
            
         // STEP 1: declare & initialize variables before ODE time-steps-------------
         //         initialize vector sizes and set all values to 0
        
         array_dims = (c_num_grps * 1 * (c_num_boxes+1)); // (3D matrix: num_grps X 1 X (num_boxes+1))
         c2_biomass_boundary_t.resize(array_dims); // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
         c_biomass_MigratorBoundary_t.resize(array_dims); // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
        
         c_external_driver_t.resize(1 * c_num_drivers * (c_num_boxes+1)); // (3D matrix: 1 X num_drivers X (num_boxes+1))
        
         c_externalForcing_grps_t.resize(1 * c_num_externalForcing_grps * c_num_boxes); // NEW!!! (3D matrix: 1 X num_externalForcing_grps X num_boxes)
                
         c_externalForcing_t.resize(1 * c_num_grps * c_num_boxes); // NEW!!! (3D matrix: 1 X num_grps X num_boxes)
                
         array_dims = (1 * c_num_grps * c_num_boxes);         // (3D matrix: 1 X num_grps X num_boxes)
         c_ConsumptionBudget_feces_t.resize(array_dims);      // (3D matrix: 1 X num_grps X num_boxes)
         c_ConsumptionBudget_metabolism_t.resize(array_dims); // (3D matrix: 1 X num_grps X num_boxes)
         c_ConsumptionBudget_eggs_t.resize(array_dims);       // (3D matrix: 1 X num_grps X num_boxes)
         c_ConsumptionBudget_predation_t.resize(array_dims);  // (3D matrix: 1 X num_grps X num_boxes)
         c_ConsumptionBudget_senescence_t.resize(array_dims); // (3D matrix: 1 X num_grps X num_boxes)
         c_ConsumptionBudget_ba_t.resize(array_dims);         // (3D matrix: 1 X num_grps X num_boxes)
         c_ConsumptionBudget_em_t.resize(array_dims);         // (3D matrix: 1 X num_grps X num_boxes)
         c_qb_t.resize(array_dims);                           // (3D matrix: 1 X num_grps X num_boxes)
         c_pb_t.resize(array_dims);                           // (3D matrix: 1 X num_grps X num_boxes)
         c_ThorntonLessemScaler_t.resize(array_dims);         // (3D matrix: 1 X num_grps X num_boxes)
         consumption_IN.resize(array_dims);                   // (3D matrix: 1 X num_grps X num_boxes)
         predation_OUT.resize(array_dims);                    // (3D matrix: 1 X num_grps X num_boxes)
         dy.resize(array_dims);                               // (3D matrix: 1 X num_grps X num_boxes) // QQQ this will get deleted when using ECOTRANode
//         dxdt.resize(array_dims); // (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes) // QQQ only needed for debugging

         flux_domain_import_driver_total_t.resize(array_dims); // (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
       
         array_dims = (c_num_grps * c_num_grps * c_num_boxes); // (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
         c_FunctionalResponse.resize(array_dims); // (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
         Q_cp_t.resize(array_dims); // (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
                
//         std::copy(std::begin(c_production_initial), std::end(c_production_initial), c_ProductionRates_t.begin()); // QQQ initialize c_ProductionRates_t


        // step 1b: pre-initialize flux variables of PhysicalFlux class ------------
        //          NOTE: pre-initialize as all 0 before the ODE solver
        //          NOTE: allows fast re-initialization to 0 at each time-step
        //          QQQ Do I really need to reset to 0 at each time-step?
        array_dims = (c_num_grps * (c_num_boxes+1)); // (2D matrix: num_grps X (num_boxes+1))
        c_preinit_flux_import_t.resize(array_dims); // declare 1 large block of memory; (2D matrix: num_grps X (num_boxes+1))
        c_preinit_flux_export_t.resize(array_dims); // declare 1 large block of memory; (2D matrix: num_grps X (num_boxes+1))
        c_preinit_flux_domain_import_t.resize(array_dims); // declare 1 large block of memory; (2D matrix: num_grps X (num_boxes+1))
        // c_preinit_flux_domain_export_t.resize(array_dims); // option for domain export values; declare 1 large block of memory; (2D matrix: num_grps X (num_boxes+1))
        
        array_dims = (1 * c_num_grps * c_num_boxes); // (3D matrix: 1 X num_grps X num_boxes DESTINY)
        c_preinit_flux_domain_import_driver_t.resize(array_dims); // declare 1 large block of memory; (3D matrix: 1 X num_grps X num_boxes DESTINY)
        c_preinit_flux_net_t.resize(array_dims); // declare 1 large block of memory; (3D matrix: 1 X num_grps X num_boxes DESTINY)
        // -------------------------------------------------------------------------

        std::cout << "Constructed class ECOTRANode" << std::endl;

    } // end constructor -----------------------------------------------------------
    

    
    // destructor ------------------------------------------------------------------
    ~ECOTRANode() {

//        std::cout << "Destroyed class ECOTRANode" << std::endl;

    } // end destructor ------------------------------------------------------------
    

    
    // method 1: m_interpolate -----------------------------------------------------
    //           interpolate time-dependent variable values @ t
    void m_interpolate(const double t_in) { // QQQ is const correct?
        t = t_in;

//        std::cout << "Running ECOTRANode.m_interpolate" << std::endl;
//        std::cout << "Class ECOTRANode: t = " << t << std::endl;

        // STEP 1: Interpolate ConsumptionBudget and physiology terms @ t ----------
        //         NOTE: this allows for defined seasonal changes in ConsumptionBudget terms
        //          3D matrix: 1 X num_grps X num_boxes)
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_grps); ++clm_index) {

                index_a     = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + clm_index; // indexing ConsumptionBudget.spline, qb.spline, & pb.spline interploaltion functions (3D matrix: 1 X num_grps X num_boxes)

                c_ConsumptionBudget_feces_t[index_a]      = c_ptr_InterpolatePhysiology[0]->spline[index_a](t);  // ConsumptionBudget_feces @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)
                c_ConsumptionBudget_metabolism_t[index_a] = c_ptr_InterpolatePhysiology[1]->spline[index_a](t);  // ConsumptionBudget_metabolism @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)
                c_ConsumptionBudget_eggs_t[index_a]       = c_ptr_InterpolatePhysiology[2]->spline[index_a](t);  // ConsumptionBudget_eggs @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)
                c_ConsumptionBudget_predation_t[index_a]  = c_ptr_InterpolatePhysiology[3]->spline[index_a](t);  // ConsumptionBudget_predation @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)
                c_ConsumptionBudget_senescence_t[index_a] = c_ptr_InterpolatePhysiology[4]->spline[index_a](t);  // ConsumptionBudget_senescence @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)
                c_ConsumptionBudget_ba_t[index_a]         = c_ptr_InterpolatePhysiology[5]->spline[index_a](t);  // ConsumptionBudget_ba @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)
                c_ConsumptionBudget_em_t[index_a]         = c_ptr_InterpolatePhysiology[6]->spline[index_a](t);  // ConsumptionBudget_em @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)
                
                c_qb_t[index_a]                           = c_ptr_InterpolatePhysiology[7]->spline[index_a](t);  // weight-specific consumption rates @ t; (fraction); (3D matrix: 1 X num_grps X num_boxes)
                c_pb_t[index_a]                           = c_ptr_InterpolatePhysiology[8]->spline[index_a](t);  // weight-specific production rates @ t; (fraction); (3D matrix: 1 X num_grps X num_boxes)
                c_ThorntonLessemScaler_t[index_a]         = c_ptr_InterpolatePhysiology[9]->spline[index_a](t);  // ThorntonLessemScaler @ t;      (fraction); (3D matrix: 1 X num_grps X num_boxes)

                // QQQ 3/9/2021 reset negative physiologies to 0, only ba & em can be negative
                std::replace_if (c_ConsumptionBudget_feces_t.begin(), c_ConsumptionBudget_feces_t.end(), [](int i){return std::signbit(i);}, 0);
                std::replace_if (c_ConsumptionBudget_metabolism_t.begin(), c_ConsumptionBudget_metabolism_t.end(), [](int i){return std::signbit(i);}, 0);
                std::replace_if (c_ConsumptionBudget_eggs_t.begin(), c_ConsumptionBudget_eggs_t.end(), [](int i){return std::signbit(i);}, 0);
                std::replace_if (c_ConsumptionBudget_predation_t.begin(), c_ConsumptionBudget_predation_t.end(), [](int i){return std::signbit(i);}, 0);
                std::replace_if (c_ConsumptionBudget_senescence_t.begin(), c_ConsumptionBudget_senescence_t.end(), [](int i){return std::signbit(i);}, 0);
                std::replace_if (c_qb_t.begin(), c_qb_t.end(), [](int i){return std::signbit(i);}, 0);
                std::replace_if (c_pb_t.begin(), c_pb_t.end(), [](int i){return std::signbit(i);}, 0);
                std::replace_if (c_ThorntonLessemScaler_t.begin(), c_ThorntonLessemScaler_t.end(), [](int i){return std::signbit(i);}, 0);
                
            } // end clm_index
        } // end layer_index
        // -------------------------------------------------------------------------

        
        // STEP 2: Interpolate physical flux rates @ t -----------------------------
        //         (m3/d); (2D matrix: num_grps X num_fluxes)
        // 3D matrix indexing is: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        for (int row_index=0; row_index<(c_num_grps); ++row_index) { // (NOTE: row_index of groups is positioned as the outer loop)
            
            // ADVECTION_flux
            for (int clm_index=0; clm_index<(c_num_fluxes[0]); ++clm_index) {
                index_a     = (0 * 1 * c_num_fluxes[0]) + (0 * c_num_fluxes[0]) + clm_index; // indexing ADVECTION_compact.spline interpolation function (2D matrix: num_t X num_fluxes_advection)
                index_b     = (0 * c_num_grps * c_num_fluxes[0]) + (row_index * c_num_fluxes[0]) + clm_index; // indexing c_FLUX_compact_t (2D matrix: num_grps X num_fluxes_advection)

                c_ptr_PhysicalFlux[0]->c_FLUX_compact_t[index_b] = c_ptr_InterpolatePhysicalFlux[0]->spline[index_a](t); // (m3/d); (2D matrix: num_grps X num_fluxes_advection); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
            } // end clm_index
            
            
            // HORIZONTALMIXING_flux
            for (int clm_index=0; clm_index<(c_num_fluxes[1]); ++clm_index) {
                index_a     = (0 * 1 * c_num_fluxes[1]) + (0 * c_num_fluxes[1]) + clm_index; // indexing HORIZONTALMIXING_compact.spline interpolation function (2D matrix: num_t X num_fluxes_HorizontalMixing)
                index_b     = (0 * c_num_grps * c_num_fluxes[1]) + (row_index * c_num_fluxes[1]) + clm_index; // indexing c_FLUX_compact_t (2D matrix: num_grps X num_fluxes_HorizontalMixing)
                
                c_ptr_PhysicalFlux[1]->c_FLUX_compact_t[index_b] = c_ptr_InterpolatePhysicalFlux[1]->spline[index_a](t); // (m3/d); (2D matrix: num_grps X num_fluxes_HorizontalMixing); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
            } // end clm_index
            
            
            // VERTICALMIXING_flux
            for (int clm_index=0; clm_index<(c_num_fluxes[2]); ++clm_index) {
                index_a     = (0 * 1 * c_num_fluxes[2]) + (0 * c_num_fluxes[2]) + clm_index; // indexing VERTICALMIXING_compact.spline interpolation function (2D matrix: num_t X num_fluxes_VerticalMixing)
                index_b     = (0 * c_num_grps * c_num_fluxes[2]) + (row_index * c_num_fluxes[2]) + clm_index; // indexing c_FLUX_compact_t (2D matrix: num_grps X num_fluxes_VerticalMixing)
                
                c_ptr_PhysicalFlux[2]->c_FLUX_compact_t[index_b] = c_ptr_InterpolatePhysicalFlux[2]->spline[index_a](t); // (m3/d); (2D matrix: num_grps X num_fluxes_VerticalMixing); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
            } // end clm_index
            
            
            // SINKING_flux
            // NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
            for (int clm_index=0; clm_index<(c_num_fluxes[3]); ++clm_index) {
                index_a     = (clm_index * 1 * c_num_grps) + (0 * c_num_grps) + row_index; // indexing SINKING_compact.spline interpolation function (3D matrix: 1 X num_grps X num_fluxes_migration)
                index_b     = (0 * c_num_grps * c_num_fluxes[3]) + (row_index * c_num_fluxes[3]) + clm_index; // indexing c_FLUX_compact_t (2D matrix: num_grps X num_fluxes_sinking)
                
               c_ptr_PhysicalFlux[3]->c_FLUX_compact_t[index_b] = c_ptr_InterpolatePhysicalFlux[3]->spline[index_a](t); // (m3/d); (2D matrix: num_grps X num_fluxes_sinking); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
            } // end clm_index
            
            
            // MIGRATION_flux
            //          NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
            for (int clm_index=0; clm_index<(c_num_fluxes[4]); ++clm_index) {
                index_a     = (clm_index * 1 * c_num_grps) + (0 * c_num_grps) + row_index; // indexing MIGRATION_compact.spline interpolation function (3D matrix: 1 X num_grps X num_fluxes_migration)
                index_b     = (0 * c_num_grps * c_num_fluxes[4]) + (row_index * c_num_fluxes[4]) + clm_index; // indexing c_FLUX_compact_t (2D matrix: num_grps X num_fluxes_migration)

                c_ptr_PhysicalFlux[4]->c_FLUX_compact_t[index_b] = c_ptr_InterpolatePhysicalFlux[4]->spline[index_a](t); // (m3/d); (2D matrix: num_grps X num_fluxes_migration); NOTE: num_fluxes is the number of flux combinations that actually exist over the whole time-series & includes domain boundary fluxes
            } // end clm_index
            
        } // end row_index (NOTE: row_index of groups is positioned as the outer loop)
        // -------------------------------------------------------------------------

        
        // STEP 3: Interpolate box dimensions @ t ----------------------------------
        //         (m3); (vector: 1 X num_boxes)
        // 3D matrix indexing is: PSUEDOARRAY[(layer_index*num_rows*num_clms) + (row_index*num_clms) + clm_index])
        for (int clm_index=0; clm_index<(c_num_boxes); ++clm_index) {
            index_a     = (0 * 1 * c_num_boxes) + (0 * c_num_boxes) + clm_index; // indexing cs_BoxVolume_t & BoxVolume.spline interpolation function (horizontal vector: 1 X num_boxes)
            PhysicalFlux::cs_BoxVolume_t[index_a] = c_ptr_InterpolatePhysicalFlux[5]->spline[index_a](t); // (m3); (vector: 1 X num_boxes)
        } // end clm_index
        // -------------------------------------------------------------------------

        
        // STEP 4: boundary biomass ------------------------------------------------
        //          NOTE: use only with NON-REFLECTIVE boundary conditions
        //          (mmoles N/m3); (3D matrix: 1 X num_grps X num_boxes+1)
        //          FFF: simplify to 2D matrix and define only for boxes with defined boundary physical fluxes
        //          QQQ gets reshaped to (3D matrix: num_grps X 1 X num_boxes+1)
        // -------------------------------------------------------------------------
        
        
        // STEP 5: interpolate external driver biomass @ t -------------------------
        //         external driver (e.g., NO3) boundary biomass conditions for each box @ t; (mmole N/m3); (3D matrix: 1 X num_drivers X (num_boxes+1))
        //         FFF: defining for all boxes now, but will try to trim down to just the boxes with fluxes (horizontal vector: 1 X num_fluxes_BoundaryImport)
        // 3D matrix indexing is: PSUEDOARRAY[(layer_index*num_rows*num_clms) + (row_index*num_clms) + clm_index])
        for (int layer_index=0; layer_index<(c_num_boxes+1); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_drivers); ++clm_index) {
                index_a     = (layer_index * 1 * c_num_drivers) + (0 * c_num_drivers) + clm_index; // indexing c_external_driver_t & external_driver.spline interpolation function (3D matrix: 1 X num_drivers X (num_boxes+1))
                c_external_driver_t[index_a]  = c_ptr_InterpolatePhysicalFlux[6]->spline[index_a](t); // external driver (e.g., NO3) boundary biomass conditions for each box @ t; (mmole N/m3); (3D matrix: 1 X num_drivers X (num_boxes+1))
            } // end for clm_index
        } // end for layer_index
        
        // QQQ 3/9/2021 reset negative physiologies to 0, only ba & em can be negative
        std::replace_if (c_external_driver_t.begin(), c_external_driver_t.end(), [](int i){return std::signbit(i);}, 0);
        // -------------------------------------------------------------------------

        
        // STEP 6: interpolate external forcing rates @ t --------------------------
        //         NEW!!! external forcing rates for each box @ t; (mmole N/m3/d); (3D matrix: 1 X num_externalForcing_grps X num_boxes)
        // 3D matrix indexing is: PSUEDOARRAY[(layer_index*num_rows*num_clms) + (row_index*num_clms) + clm_index])
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_externalForcing_grps); ++clm_index) {
                index_a     = (layer_index * 1 * c_num_externalForcing_grps) + (0 * c_num_externalForcing_grps) + clm_index; // indexing c_externalForcing_grps_t & externalForcing.spline interpolation function (3D matrix: 1 X num_externalForcing_grps X num_boxes)
                c_externalForcing_grps_t[index_a]  = c_ptr_InterpolatePhysicalFlux[7]->spline[index_a](t); // external forcing input rates (e.g., NO3) to each box @ t; (mmole N/m3); (3D matrix: 1 X num_externalForcing_grps X num_boxes)
            } // end for clm_index
        } // end for layer_index
        // -------------------------------------------------------------------------
        
        
        // STEP 7: scenario rules @ t ----------------------------------------------
        // c_Scenario_scaler_t         = QQQ; // production_input @ t; external input to model domain; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
        // -------------------------------------------------------------------------
        
        
        // STEP 8: light terms @ t -------------------------------------------------
        //          NOTE: use for Michaelis-Menten option
        // c_MLD_t                     = QQQ; // (m); (scaler)
        // c_Io_t                      = QQQ; // surface PAR light intensity @ t; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); (scaler);
        // -------------------------------------------------------------------------

    } // end method 1) m_interpolate -----------------------------------------------

    
    
//    // method qqq: m_ResetNegatives_t ---------------------------------------------------
//    //          reset negative values of c_ProductionRates_t to 0
//    //          QQQ do I still want 'const' here
//    void m_ResetNegatives_t(state_type &c_ProductionRates_t) { // QQQ void m_ResetNegatives_t(const state_type &c_ProductionRates_t)
//
//        std::replace_if (c_ProductionRates_t.begin(), c_ProductionRates_t.end(), [](int i){return std::signbit(i);}, 0);
//
//    } // end method qqq: m_ResetNegatives_t -----------------------------------------------

    
    
    
    
    
    // method 2: m_calc_Biomass_t ---------------------------------------------------
    //           calculate biomasses, nutrient concentrations, sinking flux, & boundary conditions @ t
    //           FFF Keep a running biomass time-series
    
//    void m_calc_Biomass_t(std::vector<double> ProductionRates_t_in) { // QQQ for debugging
//        c_ProductionRates_t = ProductionRates_t_in; // QQQ temp assignment from main(); c_ProductionRates_t will eventually be all internal to ECOTRANode class
    
    void m_calc_Biomass_t(const state_type &c_ProductionRates_t) {

//        std::cout << "Running ECOTRANode.m_calc_Biomass_t" << std::endl;
//        std::cout << "size c_ProductionRates_t: " << c_ProductionRates_t.size() << std::endl;

        // STEP 1: biomass @ t = q / (q/b) ---------------------------------------------
        //         3D matrix indexing: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int row_index=0; row_index<(c_num_grps); ++row_index) {
            
                index_a     = (layer_index * c_num_grps * 1) + (row_index * 1) + 0; // indexing biomass_t, biomass_boundary_t, biomass_MigratorBoundary_t & ProductionRates_t (3D matrix: num_grps X 1 X num_boxes)
                index_b     = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + row_index; // indexing qb_t (3D matrix: 1 X num_grps X num_boxes); NOTE transpose, row_index subs in for clm_index and both contain num_grps
                
                temp_biomass                            = c_ProductionRates_t[index_a] / c_qb_t[index_b]; // use qb_t value to convert rates to standing stock biomass; (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes); NOTE transpose indexing of qb_t
                
                
                // QQQ reset any negative biomass to 0
                if (std::signbit(temp_biomass)) {temp_biomass = 0;};
                
                PhysicalFlux::cs_biomass_t[index_a]     = temp_biomass; // end '+1' layer is for boundary fluxes; (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1)); NOTE access of static variable of class PhysicalFlux
                c2_biomass_boundary_t[index_a]          = temp_biomass; // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
                c_biomass_MigratorBoundary_t[index_a]   = temp_biomass; // default case is to use same boundary conditions for migration as for physical fluxes; NOTE: deactivate this line for special definition of boundary biomass for migrators (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
            
            } // end row_index
        } // end layer_index
    
        // QQQ still to code for options...
        // % biomass_PrimaryProducer_t               = biomass_t(looky_ANYPrimaryProducer, 1, :); % primary producer biomasses; (mmole N/m3); (3D matrix: num_PrimaryProducers X 1 X num_boxes); NOTE: use with Michaelis-Menten option
        // % biomass_migrator_plus1                  = biomass_migrator_t; % special definition of boundary biomass for migrators; (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes); NOTE: deactived by default
        // % biomass_migrator_plus1(:, 1, (end+1))      = 0;          % add 1 clm for boundary fluxes; (mmoles N/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: deactived by default
        // -------------------------------------------------------------------------
    
    
        // STEP 2: boundary conditions & external driver conditions ----------------
        //    OPTION 1: reflective boundary
        //              NOTE: biomass is imported into each boundary box from external environment of same biomass densities
        //              NOTE: boundary conditions are defined for each domain box
        //                    whether that box is on the edge of the domain or not
        //                    (if not, the boundary values are not used)
        // 3D matrix indexing: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        for (int layer_index=0; layer_index<(c_num_boxes+1); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_drivers); ++clm_index) {
            
                index_a         = (layer_index * 1 * c_num_drivers) + (0 * c_num_drivers) + clm_index; // indexing external_driver_t (3D matrix: 1 X num_drivers X (num_boxes+1))
                index_looky     = PhysicalFlux::cs_looky_driver[clm_index]; // row address(es) of driver group(s) (e.g., NO3); (vector: num_drivers X 1); NOTE access of static variable of class PhysicalFlux
                index_b         = (layer_index * c_num_grps * 1) + (index_looky * 1) + 0; // indexing biomass_boundary_t (3D matrix: num_grps X 1 X (num_boxes+1))
            
                c2_biomass_boundary_t[index_b] = c_external_driver_t[index_a]; // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
            
            } // end clm_index
        } // end layer_index
    
        // biomass_MigratorBoundary_t      = biomass_migrator_plus1; // NOTE: reactivate this line when using special definition of boundary biomass for migrators; (mmole N/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: deactived by default
    
    
        // %   OPTION 2: defined boundary conditions
        // %             NOTE: use for defined boundary conditions
        // %             FFF: try to trim to 2D matrix
        // biomass_boundary_t(looky_driver, 1, :)    = external_driver_t;      % paste in external (boundary) driver conditions (NO3) @ t; (mmoles/m3); (3D matrix: num_grps X 1 X num_boxes+1); NOTE: matlab automatically makes transpose of rows & clms during assignment
    
    } // end method 2) m_calc_Biomass_t ---------------------------------------------
    
    
    
    // method 3: m_calc_PhysicalFlux_t ---------------------------------------------
    //           calculate physical and migration fluxes @ t
    void m_calc_PhysicalFlux_t() {

//        std::cout << "Running ECOTRANode.m_calc_PhysicalFlux_t" << std::endl;

        // STEP 1: define variable values within PhysicalFlux class instances ------
        //          NOTE: time-dependent variables

        // ADVECTION_flux -----------
        //          NOTE: ADVECTION_flux.c_FLUX_compact_t is defined directly at the interpolation step
        std::copy(std::begin(c_preinit_flux_import_t), std::end(c_preinit_flux_import_t), c_ptr_PhysicalFlux[0]->c_flux_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_export_t), std::end(c_preinit_flux_export_t), c_ptr_PhysicalFlux[0]->c_flux_export_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_net_t), std::end(c_preinit_flux_net_t), c_ptr_PhysicalFlux[0]->c_flux_net_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
        std::copy(std::begin(c_preinit_flux_domain_import_t), std::end(c_preinit_flux_domain_import_t), c_ptr_PhysicalFlux[0]->c_flux_domain_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        // std::copy(std::begin(c_preinit_flux_domain_export_t), std::end(c_preinit_flux_domain_export_t), c_ptr_PhysicalFlux[0]->c_flux_domain_export_t.begin()); // option for domain export values; pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_domain_import_driver_t), std::end(c_preinit_flux_domain_import_driver_t), c_ptr_PhysicalFlux[0]->c_flux_domain_import_driver_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINY) NOTE: passed through ODE as ODEinput.flux_domain_import_driver_t and is not the same as external_driver_t interpolated within the ODE
        std::copy(std::begin(c2_biomass_boundary_t), std::end(c2_biomass_boundary_t), c_ptr_PhysicalFlux[0]->c_biomass_boundary_t.begin()); // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
        // --------------------------

        
        // HORIZONTALMIXING_flux ----
        //          NOTE: HORIZONTALMIXING_flux.c_FLUX_compact_t is defined directly at the interpolation step
        std::copy(std::begin(c_preinit_flux_import_t), std::end(c_preinit_flux_import_t), c_ptr_PhysicalFlux[1]->c_flux_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_export_t), std::end(c_preinit_flux_export_t), c_ptr_PhysicalFlux[1]->c_flux_export_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_net_t), std::end(c_preinit_flux_net_t), c_ptr_PhysicalFlux[1]->c_flux_net_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
        std::copy(std::begin(c_preinit_flux_domain_import_t), std::end(c_preinit_flux_domain_import_t), c_ptr_PhysicalFlux[1]->c_flux_domain_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        // std::copy(std::begin(c_preinit_flux_domain_export_t), std::end(c_preinit_flux_domain_export_t), c_ptr_PhysicalFlux[1]->c_flux_domain_export_t.begin()); // option for domain export values; pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_domain_import_driver_t), std::end(c_preinit_flux_domain_import_driver_t), c_ptr_PhysicalFlux[1]->c_flux_domain_import_driver_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINY) NOTE: passed through ODE as ODEinput.flux_domain_import_driver_t and is not the same as external_driver_t interpolated within the ODE
        std::copy(std::begin(c2_biomass_boundary_t), std::end(c2_biomass_boundary_t), c_ptr_PhysicalFlux[1]->c_biomass_boundary_t.begin()); // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
        // --------------------------
        
        
        // VERTICALMIXING_flux ------
        //          NOTE: VERTICALMIXING_flux.c_FLUX_compact_t is defined directly at the interpolation step
        std::copy(std::begin(c_preinit_flux_import_t), std::end(c_preinit_flux_import_t), c_ptr_PhysicalFlux[2]->c_flux_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_export_t), std::end(c_preinit_flux_export_t), c_ptr_PhysicalFlux[2]->c_flux_export_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_net_t), std::end(c_preinit_flux_net_t), c_ptr_PhysicalFlux[2]->c_flux_net_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
        std::copy(std::begin(c_preinit_flux_domain_import_t), std::end(c_preinit_flux_domain_import_t), c_ptr_PhysicalFlux[2]->c_flux_domain_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        // std::copy(std::begin(c_preinit_flux_domain_export_t), std::end(c_preinit_flux_domain_export_t), c_ptr_PhysicalFlux[2]->c_flux_domain_export_t.begin()); // option for domain export values; pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_domain_import_driver_t), std::end(c_preinit_flux_domain_import_driver_t), c_ptr_PhysicalFlux[2]->c_flux_domain_import_driver_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINY) NOTE: passed through ODE as ODEinput.flux_domain_import_driver_t and is not the same as external_driver_t interpolated within the ODE
        std::copy(std::begin(c2_biomass_boundary_t), std::end(c2_biomass_boundary_t), c_ptr_PhysicalFlux[2]->c_biomass_boundary_t.begin()); // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1)); QQQ do I even need boundary conditions for vertical mixing?
        // --------------------------
        
        
        // SINKING_flux -------------
        //          NOTE: SINKING_flux.c_FLUX_compact_t is defined directly at the interpolation step
        //          NOTE: there are NO boundary fluxes for sinking, c_biomass_boundary_t already defined as all 0 in constructor
        std::copy(std::begin(c_preinit_flux_import_t), std::end(c_preinit_flux_import_t), c_ptr_PhysicalFlux[3]->c_flux_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_export_t), std::end(c_preinit_flux_export_t), c_ptr_PhysicalFlux[3]->c_flux_export_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_net_t), std::end(c_preinit_flux_net_t), c_ptr_PhysicalFlux[3]->c_flux_net_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
        std::copy(std::begin(c_preinit_flux_domain_import_t), std::end(c_preinit_flux_domain_import_t), c_ptr_PhysicalFlux[3]->c_flux_domain_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        // std::copy(std::begin(c_preinit_flux_domain_export_t), std::end(c_preinit_flux_domain_export_t), c_ptr_PhysicalFlux[3]->c_flux_domain_export_t.begin()); // option for domain export values; pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_domain_import_driver_t), std::end(c_preinit_flux_domain_import_driver_t), c_ptr_PhysicalFlux[3]->c_flux_domain_import_driver_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINY) NOTE: passed through ODE as ODEinput.flux_domain_import_driver_t and is not the same as external_driver_t interpolated within the ODE
        // --------------------------
        
        
        // MIGRATION_flux -----------
        //          NOTE: MIGRATION_flux.c_FLUX_compact_t is defined directly at the interpolation step
        //          NOTE: c_biomass_boundary_t is defined as biomass_MigratorBoundary_t for MIGRATION_flux
        std::copy(std::begin(c_preinit_flux_import_t), std::end(c_preinit_flux_import_t), c_ptr_PhysicalFlux[4]->c_flux_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_export_t), std::end(c_preinit_flux_export_t), c_ptr_PhysicalFlux[4]->c_flux_export_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_net_t), std::end(c_preinit_flux_net_t), c_ptr_PhysicalFlux[4]->c_flux_net_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
        std::copy(std::begin(c_preinit_flux_domain_import_t), std::end(c_preinit_flux_domain_import_t), c_ptr_PhysicalFlux[4]->c_flux_domain_import_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        // std::copy(std::begin(c_preinit_flux_domain_export_t), std::end(c_preinit_flux_domain_export_t), c_ptr_PhysicalFlux[4]->c_flux_domain_export_t.begin()); // option for domain export values; pre-initialized as 0; (mmole N/m3/d); (2D matrix: num_grps X (num_boxes+1))
        std::copy(std::begin(c_preinit_flux_domain_import_driver_t), std::end(c_preinit_flux_domain_import_driver_t), c_ptr_PhysicalFlux[4]->c_flux_domain_import_driver_t.begin()); // pre-initialized as 0; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes DESTINY) NOTE: passed through ODE as ODEinput.flux_domain_import_driver_t and is not the same as external_driver_t interpolated within the ODE
        std::copy(std::begin(c_biomass_MigratorBoundary_t), std::end(c_biomass_MigratorBoundary_t), c_ptr_PhysicalFlux[4]->c_biomass_boundary_t.begin()); // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
        // -------------------------------------------------------------------------
        

        // STEP 2: calculate physical fluxes at t ----------------------------------
        c_ptr_PhysicalFlux[0]->m_PhysicalFlux_intraODE(); // ADVECTION_flux
        c_ptr_PhysicalFlux[1]->m_PhysicalFlux_intraODE(); // HORIZONTALMIXING_flux
        c_ptr_PhysicalFlux[2]->m_PhysicalFlux_intraODE(); // VERTICALMIXING_flux
        c_ptr_PhysicalFlux[3]->m_PhysicalFlux_intraODE(); // SINKING_flux
        c_ptr_PhysicalFlux[4]->m_PhysicalFlux_intraODE(); // MIGRATION_flux
        // -------------------------------------------------------------------------
        
        
        // STEP 3: finalize external forcing rates to each box at t ----------------
        //         NEW!!!
        // 3D matrix indexing: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_externalForcing_grps); ++clm_index) {
                
                index_a         = (layer_index * 1 * c_num_externalForcing_grps) + (0 * c_num_externalForcing_grps) + clm_index; // indexing externalForcing_grps_t (3D matrix: 1 X num_externalForcing_grps X num_boxes)
                index_looky     = PhysicalFlux::cs_looky_externalForcing[clm_index]; // clm address(es) of externally forced group(s) (e.g., NO3, juv salmon); (vector: num_externalForcing_grps X 1); NOTE access of static variable of class PhysicalFlux
                index_b         = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + index_looky; // indexing c_externalForcing_t (3D matrix: 1 X num_grps X num_boxes)
                
                c_externalForcing_t[index_b] = c_externalForcing_grps_t[index_a]; // (mmoles N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
                
            } // end clm_index
        } // end layer_index
        // -------------------------------------------------------------------------

    } // end method 3) m_calc_PhysicalFlux_t ---------------------------------------
    
    
    
    // method 4: m_calc_EnergyBudget_t ---------------------------------------------
    //           Adjust EnergyBudget to accomodate changes in ConsumptionBudget---
    //           Changes in ConsumptionBudget are due to seasonal changes in physiology, migration, etc.
    //           SSS: deactivate this step if ConsumptionBudget does not change over time
    //           NOTE: Box types are already accounted for
    //           NOTE: Make no changes to senescence due to sinking. Senescence
    //                 directs biomass transfer to detritus (which is subject to
    //                 its own sinking rate). Sinking is an additional loss term
    //                 handled as any other physical flux in dy calculation.
    //           ConsumptionBudget:
    //                       1) feces
    //                       2) metabolism
    //                       3) eggs (reproduction)
    //                       4) predation
    //                       5) senescence
    //                       6) ba (biomass accumulation)
    //                       7) em (emigration); NOTE: negative for immigration
    void m_calc_EnergyBudget_t() {

        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_grps); ++clm_index) {
                
                index_ConsumptionBudget     = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + clm_index; // indexing ConsumptionBudget_term_t (3D matrix: 1 X num_grps X num_boxes)
                
                // step 7a: ConsumptionBudget @ t: metabolism into EnergyBudget ------------
                for (int row_index=0; row_index<(c_num_nutrients); ++row_index) {
                    
                    index_fate          = (layer_index * c_num_nutrients * c_num_grps) + (row_index * c_num_grps) + clm_index; // index into fate_metabolism (3D matrix: num_nutrients X num_grps X num_boxes)
                    index_looky         = c_looky_nutrients[row_index];                                                    // QQQ what if num_nutrients=0???
                    index_EnergyBudget  = (layer_index * c_num_grps * c_num_grps) + (index_looky * c_num_grps) + clm_index;    // index into EnergyBudget (3D matrix: num_grps X num_grps X num_boxes)
                    
                    c_EnergyBudget_t[index_EnergyBudget]    = c_ConsumptionBudget_metabolism_t[index_ConsumptionBudget] * c_fate_metabolism[index_fate]; // excretion of NH4 and nitrification of NH4-->>NO3; (fraction); (3D matrix: num_nutrients X num_grps X num_boxes)
                    
                } // end for row_index
                // -------------------------------------------------------------------------
                
                
                // step 7b: ConsumptionBudget @ t: egg production into EnergyBudget --------
                for (int row_index=0; row_index<(c_num_eggs); ++row_index) {
                    
                    index_fate          = (layer_index * c_num_eggs * c_num_grps) + (row_index * c_num_grps) + clm_index; // index into fate_eggs (3D matrix: num_eggs X num_grps X num_boxes)
                    index_looky         = c_looky_eggs[row_index];                                                  // QQQ what if num_eggs=0???
                    index_EnergyBudget  = (layer_index * c_num_grps * c_num_grps) + (index_looky * c_num_grps) + clm_index; // index into EnergyBudget (3D matrix: num_grps X num_grps X num_boxes)
                    
                    c_EnergyBudget_t[index_EnergyBudget]    = c_ConsumptionBudget_eggs_t[index_ConsumptionBudget] * c_fate_eggs[index_fate]; // egg production (should work for [] eggs or for multiple eggs); (fraction); (3D matrix: num_eggs X num_grps X num_boxes)
                    
                } // end for row_index
                // -------------------------------------------------------------------------
                
                
                // step 7c: ConsumptionBudget @ t: predation & fleet catch into EnergyBudget
                for (int row_index=0; row_index<(c_num_livingANDfleets); ++row_index) {
                    
                    index_fate          = (layer_index * c_num_livingANDfleets * c_num_grps) + (row_index * c_num_grps) + clm_index; // index into fate_predation (3D matrix: num_livingANDfleets X num_grps X num_boxes)
                    index_looky         = c_looky_livingANDfleets[row_index];                                                  // QQQ what if num_livingANDfleets=0???
                    index_EnergyBudget  = (layer_index * c_num_grps * c_num_grps) + (index_looky * c_num_grps) + clm_index;    // index into EnergyBudget (3D matrix: num_grps X num_grps X num_boxes)
                    
                    c_EnergyBudget_t[index_EnergyBudget]    = c_ConsumptionBudget_predation_t[index_ConsumptionBudget] * c_fate_predation[index_fate]; // (fraction); (3D matrix: num_livingANDfleets X num_grps X num_boxes)
                    
                } // end for row_index
                // -------------------------------------------------------------------------
                
                
                // step 7d: ConsumptionBudget @ t feces & senescence into EnergyBudget -----
                for (int row_index=0; row_index<(c_num_ANYdetritus); ++row_index) {
                    
                    index_fate          = (layer_index * c_num_ANYdetritus * c_num_grps) + (row_index * c_num_grps) + clm_index; // index into fate_feces & fate_senescence (3D matrix: num_ANYdetritus X num_grps X num_boxes)
                    index_looky         = c_looky_ANYdetritus[row_index];                                                  // QQQ what if num_ANYdetritus=0???
                    index_EnergyBudget  = (layer_index * c_num_grps * c_num_grps) + (index_looky * c_num_grps) + clm_index;    // index into EnergyBudget (3D matrix: num_grps X num_grps X num_boxes)
                    
                    feces_t                             = c_ConsumptionBudget_feces_t[index_ConsumptionBudget] * c_fate_feces[index_fate];
                    senescence_t                        = c_ConsumptionBudget_senescence_t[index_ConsumptionBudget] * c_fate_senescence[index_fate];
                    c_EnergyBudget_t[index_EnergyBudget]  = feces_t + senescence_t; // detritus = feces + senescence; (fraction); (3D matrix: num_grps X num_grps X num_boxes)
                    
                } // end for row_index
                // -------------------------------------------------------------------------
                
            } // end for clm_index
        } // end for layer_index
        
    } // end method 4) m_calc_EnergyBudget_t ---------------------------------------
    
    
    
    // method 5: m_calc_CONSUMPTION_t ----------------------------------------------
    //           calculate consumption rate matrix Q_cp-----------------------------
    void m_calc_CONSUMPTION_t(const state_type &c_ProductionRates_t) {
        
        // STEP 1: arena functional responses --------------------------------------
        //          form A: Eqtn. 11 in Steele & Ruzicka 2011 ----------------------
        //              FunctionalResponseParams: producer vulnerability (m_p); (3D matrix: CONSUMERS X prey group X num_boxes) replicated across clms (= producers); NOTE: I will just use the pre-prepared 3D matrix format supplied by MATLAB and not use a vector, to allow modifications of individual trophic interctions
        //              ProductionRates_t: (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
        //              ProductionRatesC_t_repmat <--NOT USED: (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes) NOTE: VERTICAL vectors of ProductionRates_t replicated across columns; NOTE: used in MATLAB but not C++
        //              production_initial: production rates to use as initial conditions; (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
        //              productionC_initial_repmat <--NOT USED: (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: VERTICAL vectors of production_initial replicated across columns; NOTE: used in MATLAB but not C++
        //              FunctionalResponse: (matrix aligned with EnergyBudget); (unitless); (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
        // 3D matrix indexing: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_grps); ++clm_index) {          // PREY (PRODUCER) grp
                for (int row_index=0; row_index<(c_num_grps); ++row_index) {      // CONSUMER grp
                    
                    index_a     = (layer_index * c_num_grps * c_num_grps) + (row_index * c_num_grps) + clm_index; // index for FunctionalResponse & FunctionalResponseParams (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
                    index_b     = (layer_index * c_num_grps * 1) + (row_index * 1) + 0; // index for replicating clms in ProductionRates_t (=ProductionRatesC_t_repmat in MATLAB) & production_initial (= productionC_initial_repmat in MATLAB); (3D matrix: num_grps X 1 X num_boxes) NOTE: VERTICAL vectors replicated across clms
                    
                    numerator                   = ((1+c_FunctionalResponseParams[index_a]) * c_ProductionRates_t[index_b]);
                    denominator                 = (c_FunctionalResponseParams[index_a] * c_production_initial[index_b] + c_ProductionRates_t[index_b]);
                    c_FunctionalResponse[index_a] = numerator / denominator; // (matrix aligned with EnergyBudget); (unitless); (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
                    
                    if (isnan(c_FunctionalResponse[index_a])) {
                        c_FunctionalResponse[index_a] = 1; // catch NANs and set to FunctionalResponse = 1 (caused by division by 0 in 0 biomass cells); NOTE: if necessary, use isinf to +/- INF
                    }
                    
                } // end for row_index
            } // end for clm_index
        } // end for layer_index
        // -------------------------------------------------------------------------
        
        
        // STEP 2: standard calculation using ECOTRAN EnergyBudget (Acp) ----------
        //          NOTE: Q_cp is functional on its own, with or without the Michaelis-Menton option of step 7c
        //              ProductionRates_t: (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
        //              ProductionRatesP_t_repmat <--NOT USED: transpose to horizontal vectors; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: HORIZONTAL vectors of transposed ProductionRates_t replicated down rows; NOTE: used in MATLAB but not C++
        //              FunctionalResponse: (matrix aligned with EnergyBudget); (unitless); (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
        //              EnergyBudget: (unitless); (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
        //              Q_cp: (mmole N/m3/d); (3D matrix: num_grps CONSUMER X num_grps PRODUCER X num_boxes)
        // 3D matrix indexing: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_grps); ++clm_index) {          // PREY (PRODUCER) grp
                for (int row_index=0; row_index<(c_num_grps); ++row_index) {      // CONSUMER grp
                    
                    index_a         = (layer_index * c_num_grps * c_num_grps) + (row_index * c_num_grps) + clm_index; // index for Q_cp, EnergyBudget, & FunctionalResponse (3D matrix: CONSUMER num_grps X PRODUCER num_grps X num_boxes)
                    index_b         = (layer_index * c_num_grps * 1) + (clm_index * 1) + 0; // index for ProductionRates_t (3D matrix: num_grps X 1 X num_boxes); NOTE transpose clm_index for row_index
                    index_c         = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + row_index; // index for ThorntonLessemScaler_t (3D matrix: 1 X CONSUMER num_grps X num_boxes); NOTE transpose row_index for clm_index
                    
                    Q_cp_t[index_a]   = c_EnergyBudget_t[index_a] * c_FunctionalResponse[index_a] * c_ThorntonLessemScaler_t[index_c] * c_ProductionRates_t[index_b]; // (mmole N/m3/d); (3D matrix: num_grps CONSUMER X num_grps PRODUCER X num_boxes)

                } // end for row_index
            } // end for clm_index
        } // end for layer_index
        
        // QQQ 3/9/2021 reset negative Q_cp_t values to 0
        std::replace_if (Q_cp_t.begin(), Q_cp_t.end(), [](int i){return std::signbit(i);}, 0);
        // -------------------------------------------------------------------------
        
        
        // STEP 3: calculate consumption_IN & predation_OUT for each group @ t -----
        //          NOTE: predation_OUT = ProductionRates_t' if each clm of EnergyBudget sums to 1 (when FunctionalResponse = 1)
        //          NOTE: predation_OUT includes metabolism, eggs, senescence, & feces as well as predation & fleet catch
        // 3D matrix indexing: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_grps); ++clm_index) {
                
                temp_sum_clm    = 0; // re-initialize temp_sum_clm to 0 for next producer group
                temp_sum_row    = 0; // re-initialize temp_sum_rowto 0 for next consumer group
                index_c         = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + clm_index; // index into predation_OUT (3D matrix: 1 X num_grps X num_boxes) & consumption_IN (3D matrix: num_grps X 1 X num_boxes) QQQ is this difference in row and clm orientation OK to ignore?
                
                for (int row_index=0; row_index<(c_num_grps); ++row_index) {
                    
                    index_a         = (layer_index * c_num_grps * c_num_grps) + (row_index * c_num_grps) + clm_index; // for summing down columns of Q_cp
                    index_b         = (layer_index * c_num_grps * c_num_grps) + (clm_index * c_num_grps) + row_index; // for summing across rows of Q_cp; NOTE transpose of columns and rows
                    
                    temp_sum_clm    += Q_cp_t[index_a];
                    temp_sum_row    += Q_cp_t[index_b];
                    
                } // end for row_index
                
                predation_OUT[index_c]  = temp_sum_clm; // total predation flow out from each group, nutrient pool, egg pool, or detritus pool P; (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
                consumption_IN[index_c] = temp_sum_row; // total consumption flow INTO each group, nutrient pool, egg pool, or detritus pool C; (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
                
            } // end for clm_index
        } // end for layer_index
        // -------------------------------------------------------------------------
        
        
        // STEP 4: pool pelagic NH4 & benthic NH4 in all boxes ---------------------
        // 3D matrix indexing: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
        
        // step 4a: pool pelagic NH4 & benthic NH4 in all boxes
        //          NOTE: ECOTRAN model defines separate NH4 pool fates, but in this
        //                spatially-resolved format, the NH4 pool definition should be
        //                defined by the individual spatial box being considered.
        //          NOTE: code assumes there is at least 1 pelagic NH4 pool per box
        // FFF get rid of separate benthic NH4 pools in the future
        index_looky_plgcNH4 = c_looky_plgcNH4[0]; // all bnthNH4 will be added to the FIRST plgcNH4 pool, in case of multiple plgcNH4 pools per box
        
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            
            temp_sum_bnthNH4    = 0; // re-initialize temp_sum_clm to 0 for next producer group
            
            for (int clm_index=0; clm_index<(c_num_bnthNH4); ++clm_index) {
                
                index_looky_bnthNH4     = c_looky_bnthNH4[clm_index];
                index_a                 = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + index_looky_bnthNH4; // current bnthNH4 index into consumption_IN for summing bnthNH4 pools in the current box
                
                temp_sum_bnthNH4        += consumption_IN[index_a];
                consumption_IN[index_a] = 0; // zero out the current benthic NH4 pool after addding to the total benthic NH4 pool within the current box
                
            } // end for clm_index
            
            index_b                 = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + index_looky_plgcNH4; // index of first plgcNH4 pool in consumption_IN for the current box
            consumption_IN[index_b] = consumption_IN[index_b] + temp_sum_bnthNH4; // (mmole N/m3/d); (3D matrix: 1 X num_grps X num_boxes)
            
        } // end for layer_index
        // -------------------------------------------------------------------------
        
    } // end method 5) m_calc_CONSUMPTION_t ----------------------------------------


    
    // method 6: calculate dxdt ----------------------------------------------------
    //           calculate consumption rate of local input change, dxdt, from physical fluxes, migration, sinking, consumption, & predation
    //           NOTE: dy units (mmole N/m3/d2); (3D matrix: 1 X num_grps X num_boxes)
    //           NOTE: x (ProductionRates_t) units (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
    //           NOTE: use pb_t (weight-specific production rate)
    //           NOTE: physical loss terms here correspond to John's original format
    //                     where I use the term: NetAdvection_t + NetVerticalMix_t + NetHorizontalMix_t + NetSinking_t
    //                     John used the term: PhysicalLossFraction * currentProductionRate
    //                     These are the same in all essentials except that the
    //                     new code also allows for physical GAIN as well as LOSS
    //          NOTE: C++ use of ProductionRates_t for MATLAB's ProductionRates_t_transpose (this is possible because row & clm indexing are equivalent in this case)
    // 3D matrix indexing: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
    void operator() (const state_type &c_ProductionRates_t, state_type &dxdt,  double t) {
        // &x is c_ProductionRates_t
        // &dxdt is dxdt
        // t is t
        
        // STEP 1: interpolate time-dependent variable values @ t-------------------
        //         -- seven ConsumptionBudget terms
        //         -- mass-specific growth rate, pb_t
        //         -- five physical flux rates
        //         -- sub-region box dimensions
        //         -- external driver biomass
        m_interpolate(t);
        
        
//        // STEP qqq: reset negative values of c_ProductionRates_t to 0--------------
//        m_ResetNegatives_t(c_ProductionRates_t); // QQQ

        
        
        
        // STEP 2: convert ProductionRates_t to biomasses---------------------------
        m_calc_Biomass_t(c_ProductionRates_t); // QQQ do not give ProductionRates_t as a parameter because all use of ProductionRates_t will be contained within the class
        
        
        // STEP 3: calculate physical and migration fluxes @ t----------------------
        m_calc_PhysicalFlux_t();
        
        
        // STEP 4: Adjust EnergyBudget to accomodate changes in ConsumptionBudget---
        //         Changes in ConsumptionBudget are due to seasonal changes in physiology, migration, etc.
        //         SSS: deactivate this step if ConsumptionBudget does not change over time
        m_calc_EnergyBudget_t();
        
        
        // STEP 5: calculate consumption_IN & predation_OUT for each group @ t -----
        //         NOTE: predation_OUT = ProductionRates_t' if each clm of EnergyBudget sums to 1 (when FunctionalResponse = 1)
        //         NOTE: predation_OUT includes metabolism, eggs, senescence, & feces as well as predation & fleet catch
        m_calc_CONSUMPTION_t(c_ProductionRates_t);
        
        
        // STEP 6: calculate dxdt---------------------------------------------------
        for (int layer_index=0; layer_index<(c_num_boxes); ++layer_index) {
            for (int clm_index=0; clm_index<(c_num_grps); ++clm_index) {
                
                index_a       = (layer_index * 1 * c_num_grps) + (0 * c_num_grps) + clm_index; // indexing (3D matrix: 1 X num_grps X num_boxes)
                
                dxdt[index_a] = c_pb_t[index_a] * (
                                                   (c_TransferEfficiency[index_a] * (consumption_IN[index_a] + c_externalForcing_t[index_a] + c_ptr_PhysicalFlux[0]->c_flux_domain_import_driver_t[index_a] + c_ptr_PhysicalFlux[1]->c_flux_domain_import_driver_t[index_a]))
                                                   
                                                   - predation_OUT[index_a]
                                                   
                                                   - ((c_ConsumptionBudget_em_t[index_a] - c_ConsumptionBudget_ba_t[index_a]) * c_ProductionRates_t[index_a])
                                                   
                                                   // physical fluxes
                                                   + c_ptr_PhysicalFlux[0]->c_flux_net_t[index_a] // ADVECTION_flux
                                                   + c_ptr_PhysicalFlux[1]->c_flux_net_t[index_a] // HORIZONTALMIXING_flux
                                                   + c_ptr_PhysicalFlux[2]->c_flux_net_t[index_a] // VERTICALMIXING_flux
                                                   + c_ptr_PhysicalFlux[3]->c_flux_net_t[index_a] // SINKING_flux
                                                   + c_ptr_PhysicalFlux[4]->c_flux_net_t[index_a] // MIGRATION_flux
                                                   ); // NEW!!!
                
            } // end for clm_index
        } // end for layer_index
        
    } // end method 6: m_calc_dxdt_t -----------------------------------------------


}; // end class ECOTRANode
// *********************************************************************************





// observe_StateTime structure******************************************************
//      retain system state at each solution time-step
// 3D matrix indexing is: PSUEDOARRAY[(layer_index * num_rows * num_clms) + (row_index * num_clms) + clm_index])
struct observe_StateTime {

    std::vector<state_type>& s_obs_states;
    std::vector<double>& s_obs_times;

    // constructor -----------------------------------------------------------------
    observe_StateTime(std::vector<state_type> &states, std::vector<double> &times) : s_obs_states(states), s_obs_times(times) {
    } // end constructor -----------------------------------------------------------

    void operator() (const state_type &c_ProductionRates_t, double t) {
        s_obs_states.push_back(c_ProductionRates_t);
        s_obs_times.push_back(t);
    }
}; // end struct push_back_state_and_time
// *********************************************************************************
            




// declare static class variables***************************************************
// shared variables, constant over time
unsigned int PhysicalFlux::cs_num_grps; // (scalar)
unsigned int PhysicalFlux::cs_num_boxes; // (scalar)
unsigned int PhysicalFlux::cs_num_drivers; // (scalar)
std::vector<unsigned int> PhysicalFlux::cs_looky_driver; // row address(es) of driver group(s) (e.g., NO3); (vector: num_drivers X 1)
std::vector<unsigned int> PhysicalFlux::cs_looky_externalForcing; // NEW!!! column address(es) of forced external input rates of select group(s) (e.g., NO3, juv salmon); (vector: num_externalForcing_grps X 1)

// shared variables, time-dependent; NOTE: size reserved and initialized as all 0 in main
std::vector<double> PhysicalFlux::cs_biomass_t; // (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
std::vector<double> PhysicalFlux::cs_BoxVolume_t; // (m3); (vector: 1 X num_boxes)
// *********************************************************************************




 
// DEFINE THE MEX CLASS & MexFunction***********************************************
//          NOTE: All C++ MEX functions are implemented as a class named MexFunction. This class must derive from matlab::mex::Function.
class MexFunction : public matlab::mex::Function {
    
public:
    
    // All MexFunction classes must override the function call operator, operator(), to accept two arguments of type matlab::mex::ArgumentList. These arguments contain the inputs and outputs.
    void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs) {
        
        // *************************************************************************
        // STEP 1: Test to see if the arguments are of the correct type and size. If tests fail, call the MATLAB error function.
        //        checkArguments(outputs, inputs);
        // *************************************************************************
        
        
        
        
        
        // *************************************************************************
        // STEP 2: parse input arguments********************************************
        // step 2a: move input arguments into structures ---------------------------
        matlab::data::StructArray const AddressStruct_mat     = std::move(inputs[0]); // addresses & group counts; (MATLAB structure); (this also works--> matlab::data::StructArray AddressStruct_mat(inputs[0])
        matlab::data::StructArray const TrophicStruct_mat     = std::move(inputs[1]); // budgets, fates, & transfer efficiencies; (MATLAB structure)
        matlab::data::StructArray const FuncRespStruct_mat    = std::move(inputs[2]); // functional response & Michaelis Menten parameters; (MATLAB structure)
        matlab::data::StructArray const PhysicsStruct_mat     = std::move(inputs[3]); // physical fluxes & migration fluxes; (MATLAB structure)
        
        // build vectors of structure field names
        auto AddressFields      = AddressStruct_mat.getFieldNames();
        auto TrophicFields      = TrophicStruct_mat.getFieldNames();
        auto FuncRespFields     = FuncRespStruct_mat.getFieldNames();
        auto PhysicsFields      = PhysicsStruct_mat.getFieldNames();
        std::vector<matlab::data::MATLABFieldIdentifier> AddressFieldNames(AddressFields.begin(), AddressFields.end());
        std::vector<matlab::data::MATLABFieldIdentifier> TrophicFieldNames(TrophicFields.begin(), TrophicFields.end());
        std::vector<matlab::data::MATLABFieldIdentifier> FuncRespFieldNames(FuncRespFields.begin(), FuncRespFields.end());
        std::vector<matlab::data::MATLABFieldIdentifier> PhysicsFieldNames(PhysicsFields.begin(), PhysicsFields.end());
        // -------------------------------------------------------------------------
        
        
        // step 2b: get data from structure fields ---------------------------------
        //          NOTE: type definitions for matlab::data::TypedArray are defined at https://www.mathworks.com/help/matlab/apiref/
        //          NOTE: 'uint8_t' is 'unsigned char' and allows values to 255
        //          NOTE: 'uint32_t' is 'unsigned int' and allows values to 4,294,967,295

        // group counts ----------------------
        matlab::data::TypedArray<uint8_t> mat_num_grps                      = std::move(AddressStruct_mat[0][AddressFieldNames[0]]);
        matlab::data::TypedArray<uint8_t> mat_num_nutrients                 = std::move(AddressStruct_mat[0][AddressFieldNames[1]]);
        matlab::data::TypedArray<uint8_t> mat_num_NH4                       = std::move(AddressStruct_mat[0][AddressFieldNames[2]]);
        matlab::data::TypedArray<uint8_t> mat_num_plgcNH4                   = std::move(AddressStruct_mat[0][AddressFieldNames[3]]);
        matlab::data::TypedArray<uint8_t> mat_num_bnthNH4                   = std::move(AddressStruct_mat[0][AddressFieldNames[4]]);
        matlab::data::TypedArray<uint8_t> mat_num_ANYPrimaryProd            = std::move(AddressStruct_mat[0][AddressFieldNames[5]]);
        matlab::data::TypedArray<uint8_t> mat_num_phytoplankton             = std::move(AddressStruct_mat[0][AddressFieldNames[6]]);
        matlab::data::TypedArray<uint8_t> mat_num_macroalgae                = std::move(AddressStruct_mat[0][AddressFieldNames[7]]);
        matlab::data::TypedArray<uint8_t> mat_num_eggs                      = std::move(AddressStruct_mat[0][AddressFieldNames[8]]);
        matlab::data::TypedArray<uint8_t> mat_num_ANYdetritus               = std::move(AddressStruct_mat[0][AddressFieldNames[9]]);
        matlab::data::TypedArray<uint8_t> mat_num_livingANDfleets           = std::move(AddressStruct_mat[0][AddressFieldNames[10]]);
        matlab::data::TypedArray<uint8_t> mat_num_drivers                   = std::move(AddressStruct_mat[0][AddressFieldNames[11]]);
        matlab::data::TypedArray<uint8_t> mat_num_boxes                     = std::move(AddressStruct_mat[0][AddressFieldNames[12]]);
        matlab::data::TypedArray<uint32_t> mat_num_t                        = std::move(AddressStruct_mat[0][AddressFieldNames[13]]);
        matlab::data::TypedArray<uint8_t> mat_num_externalForcing_grps      = std::move(AddressStruct_mat[0][AddressFieldNames[24]]); // NEW!!! (NOTE: taking AddressFieldNames out of order here)
        // -----------------------------------

        // group addresses -------------------
        //          NOTE: ADDRESSES ARE STILL IN BASE-1 MATLAB FORMAT)
        matlab::data::TypedArray<uint8_t> mat_looky_driver                  = std::move(AddressStruct_mat[0][AddressFieldNames[14]]);
        matlab::data::TypedArray<uint8_t> mat_looky_nutrients               = std::move(AddressStruct_mat[0][AddressFieldNames[15]]);
        matlab::data::TypedArray<uint8_t> mat_looky_NO3                     = std::move(AddressStruct_mat[0][AddressFieldNames[16]]);
        matlab::data::TypedArray<uint8_t> mat_looky_NH4                     = std::move(AddressStruct_mat[0][AddressFieldNames[17]]);
        matlab::data::TypedArray<uint8_t> mat_looky_plgcNH4                 = std::move(AddressStruct_mat[0][AddressFieldNames[18]]);
        matlab::data::TypedArray<uint8_t> mat_looky_bnthNH4                 = std::move(AddressStruct_mat[0][AddressFieldNames[19]]);
        matlab::data::TypedArray<uint8_t> mat_looky_eggs                    = std::move(AddressStruct_mat[0][AddressFieldNames[20]]); // QQQ EMPTY ARRAY [] IN TEST MODEL
        matlab::data::TypedArray<uint8_t> mat_looky_livingANDfleets         = std::move(AddressStruct_mat[0][AddressFieldNames[21]]);
        matlab::data::TypedArray<uint8_t> mat_looky_ANYdetritus             = std::move(AddressStruct_mat[0][AddressFieldNames[22]]);
        matlab::data::TypedArray<uint8_t> mat_looky_externalForcing         = std::move(AddressStruct_mat[0][AddressFieldNames[23]]); // NEW!!!
        // -----------------------------------
        
        // EnergyBudget ----------------------
        matlab::data::TypedArray<double> mat_EnergyBudget                   = std::move(TrophicStruct_mat[0][TrophicFieldNames[0]]); // (3D matrix: num_grps X num_grps X num_boxes)
        // -----------------------------------

        // ConsumptionBudget terms -----------
        matlab::data::TypedArray<double> mat_ConsumptionBudget_feces        = std::move(TrophicStruct_mat[0][TrophicFieldNames[1]]); // (3D matrix: num_t X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_ConsumptionBudget_metabolism   = std::move(TrophicStruct_mat[0][TrophicFieldNames[2]]); // (3D matrix: num_t X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_ConsumptionBudget_eggs         = std::move(TrophicStruct_mat[0][TrophicFieldNames[3]]); // (3D matrix: num_t X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_ConsumptionBudget_predation    = std::move(TrophicStruct_mat[0][TrophicFieldNames[4]]); // (3D matrix: num_t X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_ConsumptionBudget_senescence   = std::move(TrophicStruct_mat[0][TrophicFieldNames[5]]); // (3D matrix: num_t X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_ConsumptionBudget_ba           = std::move(TrophicStruct_mat[0][TrophicFieldNames[6]]); // (3D matrix: num_t X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_ConsumptionBudget_em           = std::move(TrophicStruct_mat[0][TrophicFieldNames[7]]); // (3D matrix: num_t X num_grps X num_boxes)
        // -----------------------------------
        
        // TransferEfficiency ----------------
        matlab::data::TypedArray<double> mat_TransferEfficiency             = std::move(TrophicStruct_mat[0][TrophicFieldNames[8]]); // (3D matrix: 1 X num_grps X num_boxes)
        // -----------------------------------
        
        // fates -----------------------------
        matlab::data::TypedArray<double> mat_fate_feces                     = std::move(TrophicStruct_mat[0][TrophicFieldNames[9]]); // (3D matrix: num_ANYdetritus X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_fate_metabolism                = std::move(TrophicStruct_mat[0][TrophicFieldNames[10]]); // (3D matrix: num_nutrients X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_fate_eggs                      = std::move(TrophicStruct_mat[0][TrophicFieldNames[11]]); // (3D matrix: num_eggs X num_grps X num_boxes) QQQ EMPTY
        matlab::data::TypedArray<double> mat_fate_predation                 = std::move(TrophicStruct_mat[0][TrophicFieldNames[12]]); // (3D matrix: num_livingANDfleets X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_fate_senescence                = std::move(TrophicStruct_mat[0][TrophicFieldNames[13]]); // (3D matrix: num_ANYdetritus X num_grps X num_boxes)
        // -----------------------------------
        
        // initial & boundary conditions, external driver -----------
        matlab::data::TypedArray<double> mat_production_initial             = std::move(TrophicStruct_mat[0][TrophicFieldNames[14]]); // production rates to use as initial conditions; (mmole N/m3/d); (3D matrix: num_grps X 1 X num_boxes)
        matlab::data::TypedArray<double> mat_productionC_initial_repmat     = std::move(TrophicStruct_mat[0][TrophicFieldNames[15]]); // consumption inflow rate; (mmole N/m3/d); (3D matrix: num_grps X num_grps X num_boxes); NOTE: replicated vertical vectors across columns
        matlab::data::TypedArray<double> mat_external_driver                = std::move(TrophicStruct_mat[0][TrophicFieldNames[16]]); // external boundary biomass conditions for each box; external NO3 driver; (mmole N/m3); (3D matrix: num_t X num_drivers X (num_boxes+1))
        matlab::data::TypedArray<double> mat_externalForcing                = std::move(TrophicStruct_mat[0][TrophicFieldNames[17]]); // NEW!!! external input rate(s) of select group(s) to each box; (mmole N/m3/d); (3D matrix: num_t X num_externalForcing_grps X num_boxes)
        // matlab::data::TypedArray<double> mat_biomass_boundary               = std::move(TrophicStruct_mat[0][TrophicFieldNames[QQQ]]); // NOTE: use this line only for NON-REFLECTIVE boundary option; (mmoles N/m3); (3D matrix: num_t X num_grps X num_boxes+1)
        // FFF external_driver defined for all boxes now, but will try to compact to just the boxes with fluxes (3D matrix: num_t X num_drivers X num_fluxes_BoundaryImport)
        // FFF could allow for allow multiple external_driver grps?
        // -----------------------------------
        
        // scenario rules --------------------
        //          NOTE: this could be anything, using a vertical vector for scaling NO3 levels for now
        //          NNN come back to this - find dynmic scenario code in other file
        // matlab::data::TypedArray<double> mat_Scenario_scaler                = std::move(TrophicStruct_mat[0][TrophicFieldNames[QQQ]]); // scenario time-series (scaler); (3D matrix: time X num_grps X num_boxes)
        // -----------------------------------

        
        // functional response parameters ----
        matlab::data::TypedArray<double> mat_FunctionalResponseParams       = std::move(FuncRespStruct_mat[0][FuncRespFieldNames[0]]); // producer vulnerability (m_p); (3D matrix: CONSUMERS X prey group X num_boxes) replicated across clms (= producers)
        matlab::data::TypedArray<double> mat_ThorntonLessemScaler       = std::move(FuncRespStruct_mat[0][FuncRespFieldNames[1]]); // Thornton-Lessem temperature adjustment to consumption rate; value between 0 and 1; (3D matrix: num_t X num_grps X num_boxes)
        // -----------------------------------

        // optional Michaelis-Menten terms ---
        // matlab::data::TypedArray<double> mat_MichaelisMenten_Vmax           = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // Vmax  = maximum nutrient uptake rate; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_MichaelisMenten_KNO3           = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // KNO3  = NO3 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_MichaelisMenten_KNH4           = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // KNH4  = NH4 half-saturation constant; (mmol N/m3); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_MichaelisMenten_alpha          = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // alpha = initial slope of light response curve; (m2/W/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_MichaelisMenten_psi            = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // psi   = NO3 uptake inhibition by NH4; (m3/mmole N); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_MichaelisMenten_w              = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // w    = sinking rate; (m/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_MichaelisMenten_eta            = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // eta = non-grazing mortality; (1/d); (3D matrix: num_ANYPrimaryProd X 1 X num_boxes); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_MLD                            = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // time-series of Mixed Layer Depth (for light intensity calculations, NOT advection calcs); (m); (vertical vector: length = num_t); NOTE: use for Michaelis-Menten option
        // matlab::data::TypedArray<double> mat_Io                             = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // time-series of surface PAR light intensity; daily integrated solar raditation at ocean surface averaged across 24 hours; (W m^-2 h^-1); NOTE: the daily average is also what was used by Spitz in her NPZD models; NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_Kw                             = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // Light attenuation of seawater; (scaler); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_Kp                             = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // Light attenuation of phytoplankton (Newberger et al., 2003); (m2/mmol N); NOTE: use for Michaelis-Menten function
        // matlab::data::TypedArray<double> mat_ProductionFraction_macroalgae  = FuncRespStruct_mat[0][FuncRespFieldNames[QQQ]]; // tie any macroalgae production to a fixed multiple of total phytoplankton production in each box; (3D matrix: num_macroalgae X 1 X num_boxes)
        // -----------------------------------

        
        // physics terms ---------------------
        matlab::data::TypedArray<double> mat_t_grid                                         = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[0]]); // (vertical vector: num_t x 1)
        matlab::data::TypedArray<double> mat_dt                                             = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[1]]); // (vertical vector: num_t x 1)
        matlab::data::TypedArray<double> mat_BoxVolume                                      = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[2]]); // (m3); (2D matrix: num_t X num_boxes)
        matlab::data::TypedArray<double> mat_qb                                             = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[3]]); // intrinsic (weight-specific) consumption rates; (1/d); (3D matrix: num_t X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_pb                                             = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[4]]); // intrinsic (weight-specific) consumption rates; (1/d); (3D matrix: num_t X num_grps X num_boxes)
        matlab::data::TypedArray<double> mat_RetentionScaler                                = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[5]]); // (scaler 0-1); (2D matrix: num_grps X num_boxes DESTINY)
        // -----------------------------------
        
        // advection flux --------------------
        matlab::data::TypedArray<double> mat_ADVECTION_compact                              = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[6]]); // (m3/d); (2D matrix: num_t X num_fluxes_advection)
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_advection                         = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[7]]);
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_AdvectionBoundary_import          = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[8]]); //
        // matlab::data::TypedArray<uint32_t> mat_num_fluxes_AdvectionBoundary_export          = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[9]]); // option for domain export values
        matlab::data::TypedArray<double> mat_looky_AdvectionFlux                            = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[10]]); // (2D matrix: num_fluxes_advection X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        matlab::data::TypedArray<double> mat_looky_AdvectionBoundary_import                 = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[11]]); // (2D matrix: num_fluxes_AdvectionBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        // matlab::data::TypedArray<double> mat_looky_AdvectionBoundary_export                 = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[12]]); // option for domain export values; 2D matrix: num_fluxes_AdvectionBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // -----------------------------------

        // horizontal mixing flux ------------
        matlab::data::TypedArray<double> mat_HORIZONTALMIXING_compact                       = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[13]]); // (m3/d); (2D matrix: num_t X num_fluxes_HorizontalMixing)
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_HorizontalMixing                  = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[14]]);
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_HorizontalMixingBoundary_import   = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[15]]); //
        // matlab::data::TypedArray<uint32_t> mat_num_fluxes_HorizontalMixingBoundary_export   = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[16]]); // option for domain export values
        matlab::data::TypedArray<double> mat_looky_HorizontalMixingFlux                     = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[17]]); // (2D matrix: num_fluxes_HorizontalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        matlab::data::TypedArray<double> mat_looky_HorizontalMixingBoundary_import          = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[18]]); // (2D matrix: num_fluxes_HorizontalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        // matlab::data::TypedArray<double> mat_looky_HorizontalMixingBoundary_export          = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[19]]); // option for domain export values; (2D matrix: num_fluxes_HorizontalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // -----------------------------------

        // vertical mixing flux --------------
        matlab::data::TypedArray<double> mat_VERTICALMIXING_compact                         = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[20]]); // (m3/d); (2D matrix: num_t X num_fluxes_VerticalMixing)
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_VerticalMixing                    = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[21]]);
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_VerticalMixingBoundary_import     = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[22]]); //
        // matlab::data::TypedArray<uint32_t> mat_num_fluxes_VerticalMixingBoundary_export     = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[23]]); // option for domain export values
        matlab::data::TypedArray<double> mat_looky_VerticalMixingFlux                       = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[24]]); // (2D matrix: num_fluxes_VerticalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        matlab::data::TypedArray<double> mat_looky_VerticalMixingBoundary_import            = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[25]]); // (2D matrix: num_fluxes_VerticalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        // matlab::data::TypedArray<double> mat_looky_VerticalMixingBoundary_export            = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[26]]); // option for domain export values; (2D matrix: num_fluxes_VerticalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // -----------------------------------

        // sinking flux ----------------------
        matlab::data::TypedArray<double> mat_SINKING_compact                                = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[27]]); // (m3/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_sinking                           = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[28]]);
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_SinkingBoundary_import            = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[29]]); //
        // matlab::data::TypedArray<uint32_t> mat_num_fluxes_SinkingBoundary_export            = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[30]]); // option for domain export values
        matlab::data::TypedArray<double> mat_looky_SinkingFlux                              = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[31]]); // (2D matrix: num_fluxes_sinking X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        matlab::data::TypedArray<double> mat_looky_SinkingBoundary_import                   = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[32]]); // (2D matrix: num_fluxes_SinkingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        // matlab::data::TypedArray<double> mat_looky_SinkingBoundary_export                   = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[33]]); // option for domain export values; (2D matrix: num_fluxes_SinkingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // -----------------------------------
        
        // migration flux --------------------
        matlab::data::TypedArray<double> mat_MIGRATION_compact                              = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[34]]); // (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_migration                         = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[35]]);
        matlab::data::TypedArray<uint32_t> mat_num_fluxes_MigrationBoundary_import          = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[36]]); //
        // matlab::data::TypedArray<uint32_t> mat_num_fluxes_MigrationBoundary_export          = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[37]]); // option for domain export values
        matlab::data::TypedArray<double> mat_looky_MigrationFlux                            = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[38]]); // (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        matlab::data::TypedArray<double> mat_looky_MigrationBoundary_import                 = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[39]]); // (2D matrix: num_fluxes_MigrationBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        // matlab::data::TypedArray<double> mat_looky_MigrationBoundary_export                 = std::move(PhysicsStruct_mat[0][PhysicsFieldNames[40]]); // option for domain export values; (2D matrix: num_fluxes_MigrationBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // -----------------------------------

        
        // step 2c: move from matlab::data::TypedArray container to native C++ data types
        //          NOTE: must move vectors element-by-element
        // std::cout << "move from matlab::data::TypedArray container to native C++ data types" << std::endl;
        
        std::size_t array_dims;
        
        // numbers ---------------------------
        unsigned int num_grps                   = std::move(mat_num_grps[0]); // NOTE: the [0] index is needed to move from matlab::data::TypedArray
        unsigned char num_nutrients             = std::move(mat_num_nutrients[0]);
        unsigned char num_NH4                   = std::move(mat_num_NH4[0]);        // FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
        unsigned char num_plgcNH4               = std::move(mat_num_plgcNH4[0]);    // FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
        unsigned char num_bnthNH4               = std::move(mat_num_bnthNH4[0]);    // FFF replace mat_num_plgcNH4 & mat_num_bnthNH4 with simple mat_num_NH4
        unsigned char num_ANYPrimaryProd        = std::move(mat_num_ANYPrimaryProd[0]);
        unsigned char num_phytoplankton         = std::move(mat_num_phytoplankton[0]);
        unsigned char num_macroalgae            = std::move(mat_num_macroalgae[0]);
        unsigned char num_eggs                  = std::move(mat_num_eggs[0]);
        unsigned char num_ANYdetritus           = std::move(mat_num_ANYdetritus[0]);
        unsigned char num_livingANDfleets       = std::move(mat_num_livingANDfleets[0]);
        unsigned int num_drivers                = std::move(mat_num_drivers[0]);
        unsigned int num_externalForcing_grps   = std::move(mat_num_externalForcing_grps[0]); // NEW!!!
        unsigned int num_boxes                  = std::move(mat_num_boxes[0]);
        unsigned int num_t                      = std::move(mat_num_t[0]);
        double dt                               = std::move(mat_dt[0]);
        
        // put all num_fluxes_QQQ into a single num_fluxes vector
        std::vector<unsigned int> num_fluxes(10); // (vector: num_drivers X 1)
        num_fluxes[0]                       = std::move(mat_num_fluxes_advection[0]);
        num_fluxes[1]                       = std::move(mat_num_fluxes_HorizontalMixing[0]);
        num_fluxes[2]                       = std::move(mat_num_fluxes_VerticalMixing[0]);
        num_fluxes[3]                       = std::move(mat_num_fluxes_sinking[0]);
        num_fluxes[4]                       = std::move(mat_num_fluxes_migration[0]);
        num_fluxes[5]                       = std::move(mat_num_fluxes_AdvectionBoundary_import[0]);
        num_fluxes[6]                       = std::move(mat_num_fluxes_HorizontalMixingBoundary_import[0]);
        num_fluxes[7]                       = std::move(mat_num_fluxes_VerticalMixingBoundary_import[0]);
        num_fluxes[8]                       = std::move(mat_num_fluxes_SinkingBoundary_import[0]);
        num_fluxes[9]                       = std::move(mat_num_fluxes_MigrationBoundary_import[0]);
        // num_fluxes[10]                      = std::move(mat_num_fluxes_AdvectionBoundary_export[0]); // option for domain export values
        // num_fluxes[11]                      = std::move(mat_num_fluxes_HorizontalMixingBoundary_export[0]); // option for domain export values
        // num_fluxes[12]                      = std::move(mat_num_fluxes_VerticalMixingBoundary_export[0]); // option for domain export values
        // num_fluxes[13]                      = std::move(mat_num_fluxes_SinkingBoundary_export[0]); // option for domain export values
        // num_fluxes[14]                      = std::move(mat_num_fluxes_MigrationBoundary_export[0]); // option for domain export values
        // -----------------------------------
        
        
        // address variables -----------------
        //          NOTE: subtracting 1 element-by-element to put addresses into C++ base-0 format
        std::vector<unsigned int> looky_driver(num_drivers); // declare 1 large block of memory; (vector: num_drivers X 1) QQQ what if num_drivers = 0???
        for (int i=0; i < (num_drivers); ++i) {
            looky_driver[i]     = std::move(mat_looky_driver[i]) - 1; // (vector: 1 X num_drivers)
        } // end for loop; note: must move element-by-element
        
        std::vector<unsigned char> looky_nutrients(num_nutrients); // QQQ what if num_nutrients = 0???
        for (int i=0; i < (num_nutrients); ++i) {
            looky_nutrients[i]      = std::move(mat_looky_nutrients[i]) - 1; // (vector: 1 X num_nutrients)
        } // end for loop
        
        std::vector<unsigned char> looky_plgcNH4(num_plgcNH4); // QQQ what if num_plgcNH4 = 0???; FFF replace with looky_NH4
        for (int i=0; i < (num_plgcNH4); ++i) {
            looky_plgcNH4[i]        = std::move(mat_looky_plgcNH4[i]) - 1; // (vector: 1 X num_plgcNH4)
        } // end for loop
        
        std::vector<unsigned char> looky_bnthNH4(num_bnthNH4); // QQQ what if num_bnthNH4 = 0???; FFF replace with looky_NH4
        for (int i=0; i < (num_bnthNH4); ++i) {
            looky_bnthNH4[i]       = std::move(mat_looky_bnthNH4[i]) - 1; // (vector: num_bnthNH4)
        } // end for loop
        
        std::vector<unsigned char> looky_eggs(num_eggs); // QQQ what if num_eggs = 0???
        for (int i=0; i < (num_eggs); ++i) {
            looky_eggs[i]       = std::move(mat_looky_eggs[i]) - 1; // (vector: num_eggs)
        } // end for loop
        
        std::vector<unsigned char> looky_livingANDfleets(num_livingANDfleets); // QQQ what if num_livingANDfleets = 0???
        for (int i=0; i < (num_livingANDfleets); ++i) {
            looky_livingANDfleets[i]       = std::move(mat_looky_livingANDfleets[i]) - 1; // (vector: num_livingANDfleets)
        } // end for loop
        
        std::vector<unsigned char> looky_ANYdetritus(num_ANYdetritus); // QQQ what if num_ANYdetritus = 0???
        for (int i=0; i < (num_ANYdetritus); ++i) {
            looky_ANYdetritus[i]       = std::move(mat_looky_ANYdetritus[i]) - 1; // (vector: num_ANYdetritus)
        } // end for loop
        
        std::vector<unsigned int> looky_externalForcing(num_externalForcing_grps); // NEW!!! declare 1 large block of memory; (vector: num_externalForcing_grps X 1) QQQ what if num_externalForcing_grps = 0???
        for (int i=0; i < (num_externalForcing_grps); ++i) {
            looky_externalForcing[i]     = std::move(mat_looky_externalForcing[i]) - 1; // (vector: 1 X num_externalForcing_grps)
        } // end for loop; note: must move element-by-element
        // -----------------------------------

        // -----------------------------------
        std::vector<double> t_grid(num_t); // declare 1 large block of memory; (vector: 1 X num_t)
        std::move(std::begin(mat_t_grid), std::end(mat_t_grid), t_grid.begin());
        
        std::vector<double> EnergyBudget(num_grps * num_grps * num_boxes); // declare 1 large block of memory; (3D matrix: num_grps X num_grps X num_boxes)
        std::move(std::begin(mat_EnergyBudget), std::end(mat_EnergyBudget), EnergyBudget.begin());
        
        std::vector<double> FunctionalResponseParams(num_grps * num_grps * num_boxes); // declare 1 large block of memory; (3D matrix: num_grps X num_grps X num_boxes)
        std::move(std::begin(mat_FunctionalResponseParams), std::end(mat_FunctionalResponseParams), FunctionalResponseParams.begin());

        std::vector<double> TransferEfficiency(1 * num_grps * num_boxes); // declare 1 large block of memory; (3D matrix: 1 X num_grps X num_boxes)
        std::move(std::begin(mat_TransferEfficiency), std::end(mat_TransferEfficiency), TransferEfficiency.begin());

        std::vector<double> production_initial(num_grps * 1 * num_boxes); // declare 1 large block of memory; (3D matrix: num_grps X 1 X num_boxes)
        std::move(std::begin(mat_production_initial), std::end(mat_production_initial), production_initial.begin());
        
        std::vector<double> RetentionScaler(num_grps * num_boxes); // declare 1 large block of memory; (2D matrix: num_grps X num_boxes DESTINY)
        std::move(std::begin(mat_RetentionScaler), std::end(mat_RetentionScaler), RetentionScaler.begin());
        // -----------------------------------
        
        // fates -----------------------------
        std::vector<double> fate_feces(num_ANYdetritus * num_grps * num_boxes); // declare 1 large block of memory; (3D matrix: num_ANYdetritus X num_grps X num_boxes)
        std::move(std::begin(mat_fate_feces), std::end(mat_fate_feces), fate_feces.begin());

        std::vector<double> fate_metabolism(num_nutrients * num_grps * num_boxes); // declare 1 large block of memory; (3D matrix: num_nutrients X num_grps X num_boxes)
        std::move(std::begin(mat_fate_metabolism), std::end(mat_fate_metabolism), fate_metabolism.begin());

        std::vector<double> fate_eggs(num_eggs * num_grps * num_boxes); // declare 1 large block of memory; (3D matrix: num_eggs X num_grps X num_boxes)
        std::move(std::begin(mat_fate_eggs), std::end(mat_fate_eggs), fate_eggs.begin());

        std::vector<double> fate_predation(num_livingANDfleets * num_grps * num_boxes); // declare 1 large block of memory; (3D matrix: num_livingANDfleets X num_grps X num_boxes)
        std::move(std::begin(mat_fate_predation), std::end(mat_fate_predation), fate_predation.begin());

        std::vector<double> fate_senescence(num_ANYdetritus * num_grps * num_boxes); // declare 1 large block of memory; (3D matrix: num_ANYdetritus X num_grps X num_boxes)
        std::move(std::begin(mat_fate_senescence), std::end(mat_fate_senescence), fate_senescence.begin());
        // -----------------------------------
        
        // advection flux info ---------------
        array_dims = (num_fluxes[0] * 3);
        std::vector<unsigned int> looky_AdvectionFlux(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_advection X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)])
        for (int i=0; i < (array_dims); ++i) {
            looky_AdvectionFlux[i]  = std::move(mat_looky_AdvectionFlux[i]) - 1; // (2D matrix: num_fluxes_advection X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        } // end for loop
        
        array_dims = (num_fluxes[5] * 3);
        std::vector<unsigned int> looky_AdvectionBoundary_import(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        for (int i=0; i < (array_dims); ++i) {
            looky_AdvectionBoundary_import[i]   = std::move(mat_looky_AdvectionBoundary_import[i]) - 1; // (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        } // end for loop

         // option for domain export values
         // array_dims = (num_fluxes[10]* 3);
         // std::vector<unsigned int> looky_AdvectionBoundary_export(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
         // for (int i=0; i < (array_dims); ++i) {
         //     looky_AdvectionBoundary_export[i]  = std::move(mat_looky_AdvectionBoundary_export[i]) - 1; // (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
         // } // end for loop
        // -----------------------------------

        // HorizontalMixing flux info --------
        array_dims = (num_fluxes[1] * 3);
        std::vector<unsigned int> looky_HorizontalMixingFlux(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_HorizontalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)])
        for (int i=0; i < (array_dims); ++i) {
            looky_HorizontalMixingFlux[i]   = std::move(mat_looky_HorizontalMixingFlux[i]) - 1; // (2D matrix: num_fluxes_HorizontalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        } // end for loop
        
        array_dims = (num_fluxes[6] * 3);
        std::vector<unsigned int> looky_HorizontalMixingBoundary_import(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_HorizontalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        for (int i=0; i < (array_dims); ++i) {
            looky_HorizontalMixingBoundary_import[i]    = std::move(mat_looky_HorizontalMixingBoundary_import[i]) - 1; // (2D matrix: num_fluxes_HorizontalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        } // end for loop
        
        // option for domain export values
        // array_dims = (num_fluxes[11] * 3);
        // std::vector<unsigned int> looky_HorizontalMixingBoundary_export(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_HorizontalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // for (int i=0; i < (array_dims); ++i) {
        //     looky_HorizontalMixingBoundary_export[i]    = std::move(mat_looky_HorizontalMixingBoundary_export[i]) - 1; // (2D matrix: num_fluxes_HorizontalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // } // end for loop
        // -----------------------------------

        // VerticalMixing flux info ----------
        array_dims = (num_fluxes[2] * 3);
        std::vector<unsigned int> looky_VerticalMixingFlux(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_VerticalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)])
        for (int i=0; i < (array_dims); ++i) {
            looky_VerticalMixingFlux[i]   = std::move(mat_looky_VerticalMixingFlux[i]) - 1; // (2D matrix: num_fluxes_VerticalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        } // end for loop
        
        array_dims = (num_fluxes[7] * 3);
        std::vector<unsigned int> looky_VerticalMixingBoundary_import(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_VerticalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        for (int i=0; i < (array_dims); ++i) {
            looky_VerticalMixingBoundary_import[i]    = std::move(mat_looky_VerticalMixingBoundary_import[i]) - 1; // (2D matrix: num_fluxes_VerticalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        } // end for loop
        
        // option for domain export values
        // array_dims = (num_fluxes[12] * 3);
        // std::vector<unsigned int> looky_VerticalMixingBoundary_export(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_VerticalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // for (int i=0; i < (array_dims); ++i) {
        //     looky_VerticalMixingBoundary_export[i]    = std::move(mat_looky_VerticalMixingBoundary_export[i]) - 1; // (2D matrix: num_fluxes_VerticalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // } // end for loop
        // -----------------------------------

        // sinking flux info -----------------
        array_dims = (num_fluxes[3] * 3);
        std::vector<unsigned int> looky_SinkingFlux(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_sinking X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)])
        for (int i=0; i < (array_dims); ++i) {
            looky_SinkingFlux[i]   = std::move(mat_looky_SinkingFlux[i]) - 1; // (2D matrix: num_fluxes_sinking X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        } // end for loop
        
        array_dims = (num_fluxes[8] * 3);
        std::vector<unsigned int> looky_SinkingBoundary_import(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_SinkingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        for (int i=0; i < (array_dims); ++i) {
            looky_SinkingBoundary_import[i]    = std::move(mat_looky_SinkingBoundary_import[i]) - 1; // (2D matrix: num_fluxes_SinkingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        } // end for loop
        
        // option for domain export values
        // array_dims = (num_fluxes[13] * 3);
        // std::vector<unsigned int> looky_SinkingBoundary_export(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_SinkingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // for (int i=0; i < (array_dims); ++i) {
        //    looky_SinkingBoundary_export[i]    = std::move(mat_looky_SinkingBoundary_export[i]) - 1; // (2D matrix: num_fluxes_SinkingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // } // end for loop
        // -----------------------------------

        // migration flux info ---------------
        array_dims = (num_fluxes[4] * 3);
        std::vector<unsigned int> looky_MigrationFlux(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)])
        for (int i=0; i < (array_dims); ++i) {
            looky_MigrationFlux[i]   = std::move(mat_looky_MigrationFlux[i]) - 1; // (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: values constant independent of t
        } // end for loop
        
        array_dims = (num_fluxes[9] * 3);
        std::vector<unsigned int> looky_MigrationBoundary_import(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_MigrationBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        for (int i=0; i < (array_dims); ++i) {
            looky_MigrationBoundary_import[i]    = std::move(mat_looky_MigrationBoundary_import[i]) - 1; // (2D matrix: num_fluxes_MigrationBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)])
        } // end for loop
        
        // option for domain export values
        // array_dims = (num_fluxes[14] * 3);
        // std::vector<unsigned int> looky_MigrationBoundary_export(array_dims); // declare 1 large block of memory; (2D matrix: num_fluxes_MigrationBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // for (int i=0; i < (array_dims); ++i) {
        //     looky_MigrationBoundary_export[i]    = std::move(mat_looky_MigrationBoundary_export[i]) - 1; // (2D matrix: num_fluxes_MigrationBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of export flux clms in compact_flux_t)])
        // } // end for loop
        // -----------------------------------

        // *************************************************************************

        
        
        
        
        // *************************************************************************
        // STEP 3: build interpolation models---------------------------------------
        //        std::cout << "Building interpolation models" << std::endl;
        
        // step 3a: ConsumptionBudget ----------------------------------------------
        //          (3D matrix: num_t X num_grps X num_boxes)
        array_dims = (num_t * num_grps * num_boxes);
        
        std::vector<double> ConsumptionBudget_feces_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> ConsumptionBudget_metabolism_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> ConsumptionBudget_eggs_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> ConsumptionBudget_predation_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> ConsumptionBudget_senescence_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> ConsumptionBudget_em_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> ConsumptionBudget_ba_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> qb_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> pb_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        std::vector<double> ThorntonLessemScaler_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_boxes DESTINY)
        
        std::move(std::begin(mat_ConsumptionBudget_feces), std::end(mat_ConsumptionBudget_feces), ConsumptionBudget_feces_series.begin());
        std::move(std::begin(mat_ConsumptionBudget_metabolism), std::end(mat_ConsumptionBudget_metabolism), ConsumptionBudget_metabolism_series.begin());
        std::move(std::begin(mat_ConsumptionBudget_eggs), std::end(mat_ConsumptionBudget_eggs), ConsumptionBudget_eggs_series.begin());
        std::move(std::begin(mat_ConsumptionBudget_predation), std::end(mat_ConsumptionBudget_predation), ConsumptionBudget_predation_series.begin());
        std::move(std::begin(mat_ConsumptionBudget_senescence), std::end(mat_ConsumptionBudget_senescence), ConsumptionBudget_senescence_series.begin());
        std::move(std::begin(mat_ConsumptionBudget_ba), std::end(mat_ConsumptionBudget_ba), ConsumptionBudget_ba_series.begin());
        std::move(std::begin(mat_ConsumptionBudget_em), std::end(mat_ConsumptionBudget_em), ConsumptionBudget_em_series.begin());
        std::move(std::begin(mat_qb), std::end(mat_qb), qb_series.begin()); // intrinsic (weight-specific) consumption rates; (1/d); (3D matrix: num_t X num_grps X num_boxes)
        std::move(std::begin(mat_pb), std::end(mat_pb), pb_series.begin()); // intrinsic (weight-specific) production rates; (1/d); (3D matrix: num_t X num_grps X num_boxes)
        std::move(std::begin(mat_ThorntonLessemScaler), std::end(mat_ThorntonLessemScaler), ThorntonLessemScaler_series.begin()); // (3D matrix: num_t X num_grps X num_boxes)

        SplineInterpolate ConsumptionBudget_feces(ConsumptionBudget_feces_series,           num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate ConsumptionBudget_metabolism(ConsumptionBudget_metabolism_series, num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate ConsumptionBudget_eggs(ConsumptionBudget_eggs_series,             num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate ConsumptionBudget_predation(ConsumptionBudget_predation_series,   num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate ConsumptionBudget_senescence(ConsumptionBudget_senescence_series, num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate ConsumptionBudget_ba(ConsumptionBudget_ba_series,                 num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate ConsumptionBudget_em(ConsumptionBudget_em_series,                 num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate qb(qb_series,                                                     num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate pb(pb_series,                                                     num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        SplineInterpolate ThorntonLessemScaler(ThorntonLessemScaler_series,                 num_t, 1, num_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_grps X num_boxes)
        // -------------------------------------------------------------------------
        
        
        // step 3b: physical flux rates --------------------------------------------
        // advection flux
        array_dims = (num_t * num_fluxes[0]);
        std::vector<double> ADVECTION_compact_series(array_dims); // declare 1 large block of memory; (2D matrix: num_t X num_fluxes_advection)
        std::move(std::begin(mat_ADVECTION_compact), std::end(mat_ADVECTION_compact), ADVECTION_compact_series.begin()); // (m3/d); (2D matrix: num_t X num_fluxes_advection)
        SplineInterpolate ADVECTION_compact(ADVECTION_compact_series, num_t, 1, num_fluxes[0], 1, dt); // interpolation function; (horizontal vector: 1 X num_fluxes_advection)
        
        // HorizontalMixing flux
        array_dims = (num_t * num_fluxes[1]);
        std::vector<double> HORIZONTALMIXING_compact_series(array_dims); // declare 1 large block of memory; (2D matrix: num_t X num_fluxes_HorizontalMixing)
        std::move(std::begin(mat_HORIZONTALMIXING_compact), std::end(mat_HORIZONTALMIXING_compact), HORIZONTALMIXING_compact_series.begin()); // (m3/d); (2D matrix: num_t X num_fluxes_HorizontalMixing)
        SplineInterpolate HORIZONTALMIXING_compact(HORIZONTALMIXING_compact_series, num_t, 1, num_fluxes[1], 1, dt); // interpolation function; (horizontal vector: 1 X num_fluxes_HorizontalMixing)
        
        // VerticalMixing flux
        array_dims = (num_t * num_fluxes[2]);
        std::vector<double> VERTICALMIXING_compact_series(array_dims); // declare 1 large block of memory; (2D matrix: num_t X num_fluxes_VerticalMixing)
        std::move(std::begin(mat_VERTICALMIXING_compact), std::end(mat_VERTICALMIXING_compact), VERTICALMIXING_compact_series.begin()); // (m3/d); (2D matrix: num_t X num_fluxes_VerticalMixing)
        SplineInterpolate VERTICALMIXING_compact(VERTICALMIXING_compact_series, num_t, 1, num_fluxes[2], 1, dt); // interpolation function; (horizontal vector: 1 X num_fluxes_VerticalMixing)

        // sinking flux (NOTE: 3D matrix)
        array_dims = (num_t * num_grps * num_fluxes[3]);
        std::vector<double> SINKING_compact_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_fluxes_sinking)
        std::move(std::begin(mat_SINKING_compact), std::end(mat_SINKING_compact), SINKING_compact_series.begin()); // (m3/d); (3D matrix: num_t X num_grps X num_fluxes_sinking)
        SplineInterpolate SINKING_compact(SINKING_compact_series, num_t, 1, num_grps, num_fluxes[3], dt); // interpolation function; (3D matrix: 1 X num_grps X num_fluxes_sinking)

        // migration flux (NOTE: 3D matrix)
        array_dims = (num_t * num_grps * num_fluxes[4]);
        std::vector<double> MIGRATION_compact_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_grps X num_fluxes_migration)
        std::move(std::begin(mat_MIGRATION_compact), std::end(mat_MIGRATION_compact), MIGRATION_compact_series.begin()); // (m3/d); (3D matrix: num_t X num_grps X num_fluxes_migration)
        SplineInterpolate MIGRATION_compact(MIGRATION_compact_series, num_t, 1, num_grps, num_fluxes[4], dt); // interpolation function; (3D matrix: 1 X num_grps X num_fluxes_migration)
        // -------------------------------------------------------------------------


        // step 3c: box volume -----------------------------------------------------
        array_dims = (num_t * num_boxes);
        std::vector<double> BoxVolume_series(array_dims); // declare 1 large block of memory; (2D matrix: num_t X num_boxes)
        std::move(std::begin(mat_BoxVolume), std::end(mat_BoxVolume), BoxVolume_series.begin()); // (m3); (2D matrix: num_t X num_boxes)
        SplineInterpolate BoxVolume(BoxVolume_series, num_t, 1, num_boxes, 1, dt); // interpolation function; (horizontal vector: 1 X num_boxes)
        // -------------------------------------------------------------------------

        
        // step 3d: boundary biomass -----------------------------------------------
        //          NOTE: use only with NON-REFLECTIVE boundary conditions
        //          FFF: simplify to 2D matrix and define only for boxes with defined boundary physical fluxes
        // -------------------------------------------------------------------------


        // step 3e: external driver biomass ----------------------------------------
        //          external NO3 driver; external boundary biomass conditions for each box; (mmole N/m3); (3D matrix: num_t X num_drivers X (num_boxes+1))
        array_dims = (num_t * num_drivers * (num_boxes+1));
        std::vector<double> external_driver_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_drivers X (num_boxes+1))
        std::move(std::begin(mat_external_driver), std::end(mat_external_driver), external_driver_series.begin()); // external boundary biomass conditions for each box; external NO3 driver; (mmole N/m3); (3D matrix: num_t X num_drivers X (num_boxes+1))
        SplineInterpolate external_driver(external_driver_series, num_t, 1, num_drivers, (num_boxes+1), dt); // interpolation function; (3D matrix: 1 X num_drivers X (num_boxes+1))
        // -------------------------------------------------------------------------
        
        
        // step 3f: externalForcing rates -------------------------------------------
        //          NEW!!! forced input rates to each box; (mmole N/m3/d); (3D matrix: num_t X num_externalForcing_grps X num_boxes)
        array_dims = (num_t * num_externalForcing_grps * num_boxes);
        std::vector<double> externalForcing_series(array_dims); // declare 1 large block of memory; (3D matrix: num_t X num_externalForcing_grps X num_boxes)
        std::move(std::begin(mat_externalForcing), std::end(mat_externalForcing), externalForcing_series.begin()); // forced input rates to each box; (mmole N/m3/d); (3D matrix: num_t X num_externalForcing_grps X num_boxes)
        SplineInterpolate externalForcing(externalForcing_series, num_t, 1, num_externalForcing_grps, num_boxes, dt); // interpolation function; (3D matrix: 1 X num_externalForcing_grps X num_boxes)
        // -------------------------------------------------------------------------
        
        
        // step 3g: scenario rules -------------------------------------------------
        // -------------------------------------------------------------------------
        
        
        // step 3h: light terms ----------------------------------------------------
        //          NOTE: use for Michaelis-Menten option
        // -------------------------------------------------------------------------
        
        
        // step 3i: build pointer vectors to interpolation functions ---------------
        // ConsumptionBudget
        std::vector<SplineInterpolate*> ptr_InterpolatePhysiology; // declare vector of pointers
        ptr_InterpolatePhysiology.resize(10); // initialize vector size to all 0; (vector: 10 X 1)
        ptr_InterpolatePhysiology[0] = &ConsumptionBudget_feces; // assign the addresses of each SplineInterpolate ConsumptionBudget instance
        ptr_InterpolatePhysiology[1] = &ConsumptionBudget_metabolism;
        ptr_InterpolatePhysiology[2] = &ConsumptionBudget_eggs;
        ptr_InterpolatePhysiology[3] = &ConsumptionBudget_predation;
        ptr_InterpolatePhysiology[4] = &ConsumptionBudget_senescence;
        ptr_InterpolatePhysiology[5] = &ConsumptionBudget_ba;
        ptr_InterpolatePhysiology[6] = &ConsumptionBudget_em;
        ptr_InterpolatePhysiology[7] = &qb;
        ptr_InterpolatePhysiology[8] = &pb;
        ptr_InterpolatePhysiology[9] = &ThorntonLessemScaler; // assign the address of SplineInterpolate ThorntonLessemScaler instance

        // physical fluxes & other physical terms
        std::vector<SplineInterpolate*> ptr_InterpolatePhysicalFlux; // declare vector of pointers
        ptr_InterpolatePhysicalFlux.resize(8); // initialize vector size to all 0; (vector: 8 X 1)
        ptr_InterpolatePhysicalFlux[0] = &ADVECTION_compact; // assign the addresses of each SplineInterpolate physical flux instance
        ptr_InterpolatePhysicalFlux[1] = &HORIZONTALMIXING_compact;
        ptr_InterpolatePhysicalFlux[2] = &VERTICALMIXING_compact;
        ptr_InterpolatePhysicalFlux[3] = &SINKING_compact;
        ptr_InterpolatePhysicalFlux[4] = &MIGRATION_compact;
        ptr_InterpolatePhysicalFlux[5] = &BoxVolume;
        ptr_InterpolatePhysicalFlux[6] = &external_driver;
        ptr_InterpolatePhysicalFlux[7] = &externalForcing; // NEW!!!
        // FFF add pointer to boundary biomass when using NON-REFLECTIVE boundary conditions

        // FFF scenario rules
        
        // FFF light terms when using Michaelis-Menten option
        
        // *************************************************************************
        
        
        

        
        // *************************************************************************
        // STEP 5: initialize PhysicalFlux classes----------------------------------
        // step 5a: initialize static member of PhysicalFlux class -----------------
        // initialize static (shared) vector size to all 0
        PhysicalFlux::cs_looky_driver.resize(num_drivers); // initialize vector size to all 0; row address(es) of driver group(s) (e.g., NO3); (vector: num_drivers)
        PhysicalFlux::cs_biomass_t.resize(num_grps * 1 * (num_boxes+1)); // initialize vector size to all 0; (mmoles N/m3); (3D matrix: num_grps X 1 X (num_boxes+1))
        PhysicalFlux::cs_BoxVolume_t.resize(num_boxes); // initialize vector size to all 0; (m3); (vector: 1 X num_boxes)
        PhysicalFlux::cs_looky_externalForcing.resize(num_externalForcing_grps); // NEW!!! initialize vector size to all 0; row address(es) of external forcing group(s) (e.g., NO3, juv salmon); (vector: num_externalForcing_grps)

        // define shared variable values, constant over time
        PhysicalFlux::cs_num_grps          = num_grps;
        PhysicalFlux::cs_num_boxes         = num_boxes;
        PhysicalFlux::cs_num_drivers       = num_drivers;
        std::copy(std::begin(looky_driver), std::end(looky_driver), PhysicalFlux::cs_looky_driver.begin()); // row address(es) of driver group(s) (e.g., NO3); (vector: num_drivers); NOTE use of copy not move
        std::copy(std::begin(looky_externalForcing), std::end(looky_externalForcing), PhysicalFlux::cs_looky_externalForcing.begin()); // NEW!!! row address(es) of external forcing group(s) (e.g., NO3, juv salmon); (vector: num_externalForcing_grps); NOTE use of copy not move
        // -------------------------------------------------------------------------

        
        // step 5b: instantiate ADVECTION_flux -------------------------------------
        PhysicalFlux ADVECTION_flux(num_fluxes[0], num_fluxes[5] /* , num_fluxes[10] */); // instantiate ADVECTION_flux with constructor; NOTE option for domain export values
        
        // flux-specific variables, constant over time
        std::move(std::begin(looky_AdvectionFlux), std::end(looky_AdvectionFlux), ADVECTION_flux.c_looky_flux.begin()); // (2D matrix: num_fluxes X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: vector size already defined in constructor
        std::move(std::begin(looky_AdvectionBoundary_import), std::end(looky_AdvectionBoundary_import), ADVECTION_flux.c_looky_boundary_import.begin()); // (2D matrix: num_fluxes_BoundaryImport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        // std::move(std::begin(looky_AdvectionBoundary_export), std::end(looky_AdvectionBoundary_export), ADVECTION_flux.c_looky_boundary_export.begin()); // option for domain export values; (2D matrix: num_fluxes_BoundaryExport X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        std::move(std::begin(RetentionScaler), std::end(RetentionScaler), ADVECTION_flux.c_RetentionScaler.begin()); // (2D matrix: num_grps X num_boxes DESTINY)
        // -------------------------------------------------------------------------

        
        // step 5c: instantiate HORIZONTALMIXING_flux ------------------------------
        PhysicalFlux HORIZONTALMIXING_flux(num_fluxes[1], num_fluxes[6] /* , num_fluxes_[11] */); // instantiate HORIZONTALMIXING_flux with constructor; NOTE option for domain export values
        
        // flux-specific variables, constant over time
        std::move(std::begin(looky_HorizontalMixingFlux), std::end(looky_HorizontalMixingFlux), HORIZONTALMIXING_flux.c_looky_flux.begin()); // (2D matrix: num_fluxes_HorizontalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: vector size already defined in constructor
        std::move(std::begin(looky_HorizontalMixingBoundary_import), std::end(looky_HorizontalMixingBoundary_import), HORIZONTALMIXING_flux.c_looky_boundary_import.begin()); // (2D matrix: num_fluxes_HorizontalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        // std::move(std::begin(looky_HorizontalMixingBoundary_export), std::end(looky_HorizontalMixingBoundary_export), HORIZONTALMIXING_flux.c_looky_boundary_export.begin()); // option for domain export values; (2D matrix: num_fluxes_HorizontalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        std::move(std::begin(RetentionScaler), std::end(RetentionScaler), HORIZONTALMIXING_flux.c_RetentionScaler.begin()); // (2D matrix: num_grps X num_boxes DESTINY)
        // -------------------------------------------------------------------------

        
        // step 5d: instantiate VERTICALMIXING_flux ------------------------------
        PhysicalFlux VERTICALMIXING_flux(num_fluxes[2], num_fluxes[7] /*, num_fluxes[12] */); // instantiate VERTICALMIXING_flux with constructor; NOTE option for domain export values
        
        // flux-specific variables, constant over time
        std::move(std::begin(looky_VerticalMixingFlux), std::end(looky_VerticalMixingFlux), VERTICALMIXING_flux.c_looky_flux.begin()); // (2D matrix: num_fluxes_VerticalMixing X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: vector size already defined in constructor
        std::move(std::begin(looky_VerticalMixingBoundary_import), std::end(looky_VerticalMixingBoundary_import), VERTICALMIXING_flux.c_looky_boundary_import.begin()); // (2D matrix: num_fluxes_VerticalMixingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        // std::move(std::begin(looky_VerticalMixingBoundary_export), std::end(looky_VerticalMixingBoundary_export), VERTICALMIXING_flux.c_looky_boundary_export.begin()); // option for domain export values; (2D matrix: num_fluxes_VerticalMixingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        std::move(std::begin(RetentionScaler), std::end(RetentionScaler), VERTICALMIXING_flux.c_RetentionScaler.begin()); // (2D matrix: num_grps X num_boxes DESTINY)
        // -------------------------------------------------------------------------

 
        // step 5e: instantiate SINKING_flux ---------------------------------------
        //          NOTE: SINKING_flux.c_RetentionScaler values are always 0 and already defined in contructor
        PhysicalFlux SINKING_flux(num_fluxes[3], num_fluxes[8] /* , num_fluxes[13] */); // instantiate SINKING_flux with constructor; NOTE option for domain export values
        
        // flux-specific variables, constant over time
        std::move(std::begin(looky_SinkingFlux), std::end(looky_SinkingFlux), SINKING_flux.c_looky_flux.begin()); // (2D matrix: num_fluxes_sinking X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: vector size already defined in constructor
        std::move(std::begin(looky_SinkingBoundary_import), std::end(looky_SinkingBoundary_import), SINKING_flux.c_looky_boundary_import.begin()); // (2D matrix: num_fluxes_SinkingBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        // std::move(std::begin(looky_SinkingBoundary_export), std::end(looky_SinkingBoundary_export), SINKING_flux.c_looky_boundary_export.begin()); // option for domain export values; (2D matrix: num_fluxes_SinkingBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        // -------------------------------------------------------------------------

        
        // step 5f: instantiate MIGRATION_flux -------------------------------------
        //          NOTE: MIGRATION_flux.c_RetentionScaler values are always 0 and already defined in contructor
        PhysicalFlux MIGRATION_flux(num_fluxes[4], num_fluxes[9] /* , num_fluxes[14] */); // instantiate MIGRATION_flux with constructor; NOTE option for domain export values
        
        // flux-specific variables, constant over time
        std::move(std::begin(looky_MigrationFlux), std::end(looky_MigrationFlux), MIGRATION_flux.c_looky_flux.begin()); // (2D matrix: num_fluxes_migration X 3-->[(DESTINY box) (SOURCE box) (flux address in non-compacted flux 2D matrix: DESTINY X SOURCE)]); NOTE: fluxes include all linked boxes +1 for external links; NOTE: vector size already defined in constructor
        std::move(std::begin(looky_MigrationBoundary_import), std::end(looky_MigrationBoundary_import), MIGRATION_flux.c_looky_boundary_import.begin()); // (2D matrix: num_fluxes_MigrationBoundary_import X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        // std::move(std::begin(looky_MigrationBoundary_export), std::end(looky_MigrationBoundary_export), MIGRATION_flux.c_looky_boundary_export.begin()); // option for domain export values; (2D matrix: num_fluxes_MigrationBoundary_export X 3-->[(DESTINY box) (SOURCE box) (addresses of import flux clms in compact_flux_t)]); NOTE: vector size already defined in constructor
        // -------------------------------------------------------------------------
        
        
        // step 5g: create pointers to PhysicalFlux instances ----------------------
        std::vector<PhysicalFlux*> ptr_PhysicalFlux; // declare vector of pointers of type PhysicalFlux
        ptr_PhysicalFlux.resize(5); // initialize vector size to all 0; (vector: 5 X 1)
        ptr_PhysicalFlux[0] = &ADVECTION_flux;
        ptr_PhysicalFlux[1] = &HORIZONTALMIXING_flux;
        ptr_PhysicalFlux[2] = &VERTICALMIXING_flux;
        ptr_PhysicalFlux[3] = &SINKING_flux;
        ptr_PhysicalFlux[4] = &MIGRATION_flux;
        // *************************************************************************

        
        
        
        
        
        
        
        // *************************************************************************
        // STEP 6: SOLVE THE ODE----------------------------------------------------
        // step 6a: instantiate ODE class ------------------------------------------
        ECOTRANode E2Etest(num_grps, num_nutrients, num_bnthNH4, num_eggs, num_ANYdetritus, num_livingANDfleets, num_boxes, num_drivers, num_externalForcing_grps, num_fluxes, looky_nutrients, looky_eggs, looky_livingANDfleets, looky_ANYdetritus, looky_plgcNH4, looky_bnthNH4, ptr_InterpolatePhysiology, ptr_InterpolatePhysicalFlux, ptr_PhysicalFlux, EnergyBudget, fate_feces, fate_metabolism, fate_eggs, fate_predation, fate_senescence, TransferEfficiency, FunctionalResponseParams, production_initial); // NEW!!!
        // -------------------------------------------------------------------------

        
        // step 6b: initialize the ODE ---------------------------------------------
        state_type c_ProductionRates_t; // (3D matrix: num_grps X 1 X num_boxes DESTINY) QQQ is this the place to declare ProductionRates_t? QQQ change name to ProductionRates_t
        c_ProductionRates_t.resize(num_grps * 1 * num_boxes); //  initialize ProductionRates_t as all 0; (3D matrix: num_grps X 1 X num_boxes DESTINY)
        std::copy(std::begin(production_initial), std::end(production_initial), c_ProductionRates_t.begin()); // initialize ProductionRates_t

        // define start and end times for odeint
        double t_start, t_end;
        t_start     = t_grid[0];
//        t_end       = t_grid[(num_t - 1)] + 1; // add 1 to run odeint solver through the final day; QQQ ERROR runs 1 day too long
        
//        t_end       = t_grid[(num_t)] + 1; // QQQ ERROR returns empty variables
//        t_end       = t_grid[(num_t - 1)]; // QQQ ERROR 1 time-step shorter than t_grid
        t_end       = t_grid[(num_t - 1)] + dt; // 1 time-step shorter than t_grid


        
        
        
        std::cout << "just initialized c_ProductionRates_t" << std::endl; // QQQ
//        std::cout << "total size: " << c_ProductionRates_t.size() << std::endl; // QQQ
        // -------------------------------------------------------------------------

        
        // step 6c: solve the ODE for each time-step--------------------------------
        std::cout << "STARTING ODE" << std::endl;

        // declare the time-series observer variables
        std::vector<state_type> obs_states; // rename as re_Y
        std::vector<double> obs_times;
        size_t steps;
        
        steps = boost::numeric::odeint::integrate_const(stepper, E2Etest, c_ProductionRates_t, t_start, t_end, dt, observe_StateTime(obs_states, obs_times)); // use observer to return the solution for x at each time-step; WORKS!!! 10/23/2020

        std::cout << "FINISHED ODE" << std::endl;
        // *************************************************************************
        
        
        
        
        
        // *************************************************************************
        // STEP 7: output results to MATLAB for debugging---------------------------
        std::cout << "export results to MATLAB" << std::endl;
        
        
        // step 7a: copy C++ variables to matlab variables -------------------------
        size_t num_clms = (num_grps * 1 * num_boxes);
        matlab::data::TypedArray<double> mat_obs_times      = factory.createArray<double>({steps});
        matlab::data::TypedArray<double> mat_obs_states     = factory.createArray<double>({steps, num_clms});
        
        // use loop to copy element-by-element obs_times & obs_states to mat_obs_times & mat_obs_states (THIS WORKS)
        for (size_t i(0); i < steps; ++i) {
            mat_obs_times[i]    = obs_times[i];
            for (size_t j(0); j < num_clms; ++j) {
                mat_obs_states[i][j] = obs_states[i][j];
            } // end for j (clms)
        } // end for i (rows)
        // -------------------------------------------------------------------------
        
        
        // step 7b: define 2-element struct array with two fields per struct -------
        auto structArray = factory.createStructArray({1}, {
            
            "mat_obs_times",
            "mat_obs_states",
            
        });
        // -------------------------------------------------------------------------

        
        // step 7c: assign values to each field in first struct --------------------
        structArray[0]["mat_obs_times"]             = mat_obs_times;
        structArray[0]["mat_obs_states"]            = mat_obs_states;
        // -------------------------------------------------------------------------



        // step 7d: copy structArray to MEX variable :outputs" ---------------------
        outputs[0] = structArray;
        // *************************************************************************

        
    } // end void operator()(matlab::mex::ArgumentList outputs, matlab::mex::ArgumentList inputs)
}; // end class MexFunction : public matlab::mex::Function

// end file ************************************************************************
