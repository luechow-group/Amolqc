cmake_minimum_required(VERSION 3.10)

add_subdirectory(minimizers)

add_fortran_library(src STATIC
        # Global
        global_m.f90
        machine_m.f90
        init_m.f90
        mainLoop_m.f90
        subloop_m.f90
        # E_loc
        waveFunction_m.f90
        wfData_m.f90
        jastrowDTN_m.f90
        jastrowIC_m.f90
        jastrowAniso_m.f90
        jastrowParamData_m.f90
        ecpIo_m.f90
        ecp_m.f90
        fPsi2_m.f90
        multiDet_m.f90
        multiDetParam_m.f90
        moParam_m.f90
        aoMoTask_m.f90
        eConfigs_m.f90
        linkedList_m.f90
        refUtils_m.f90
        refBase_m.f90
        refSimple_m.f90
        refStr_m.f90
        refStrPerm_m.f90
        refVal_m.f90
        refValStr_m.f90
        refPos_m.f90
        refCtr_m.f90
        refVector_m.f90
        refADT_m.f90
        refVList_m.f90
        refListVList_m.f90
        refDomain_m.f90
        references_m.f90
        findNucElecs_m.f90
        rdataupdate_m.f90
        sphericalIntGrids_m.f90
        jastrow_m.f90
        jastrowSM_m.f90
        aos_m.f90
        aosData_m.f90
        cuspOpt_m.f90
        elocData_m.f90
        eloc_m.f90
        mos_m.f90
        aoCut.f90
        aoMo_m.f90
        aoMoCut.f90
        elocTest_m.f90
        coulombDensity_m.f90
        # QMC
        qmc_m.f90
        qmcSample_m.f90
        propagator_m.f90
        randomWalker_m.f90
        rwSample_m.f90
        rwStatistics_m.f90
        reconf_m.f90
        properties_m.f90
        initialPositions_m.f90
        # Wave Function Optimization
        optimizeParams_m.f90
        wfParameters_m.f90
        optDerivsTest_m.f90
        elocAndPsiTermsBase_m.f90
        elocAndPsiTermsLin_m.f90
        elocAndPsiTermsLM_m.f90
        elocAndPsiTermsWEBFGS_m.f90
        elocAndPsiTermsENR_m.f90
        elocAndPsiTermsEBFGS_m.f90
        elocAndPsiTermsGen_m.f90
        optParamsBFGS_m.f90
        optParamsELin_m.f90
        optParamsENR_m.f90
        optParamsLBFGS_m.f90
        optParamsPOpt_m.f90
        optParamsVarmin_m.f90
        optParamsWBFGS_m.f90
        optParamsVNL2SOL_m.f90
        # Wave Function Analysis
        eigenVectAnalysis_m.f90
        psiMax_m.f90
        fGradPsi2_m.f90
        #sed_m.f90
        assign_m.f90
        findNucElecs_m.f90
        energyPart_m.f90
        hungarian_m.f90
        maxBasins_m.f90
        maxAnalysis_m.f90
        maxRawData_m.f90
        maximizeSample_m.f90
        electronDensity_m.f90
        rhoMax_m.f90
        rhoData_m.f90
        partition_m.f90
        rhoGrid_m.f90
        maximizeSampleRho_m.f90
        moMax_m.f90
        posList_m.f90
        vxc_m.f90
        # MPI
        mpiInterface_m.f90)

set_file_compile_options(aoMoTask_m.f90 ${CMAKE_Fortran_OMP_FLAG})

set_file_compile_options(init_m.f90 ${VERSION_FLAG})

set_file_compile_options(jastrowIC_m.f90 ${VECTORINFO_FLAGS})
set_file_compile_options(mos_m.f90 ${VECTORINFO_FLAGS})

if (${CHKNANUP})
    set_file_compile_options(multiDet_m.f90 "-DCHKNANUP")
    set_file_compile_options(waveFunction_m.f90 "-DCHKNANUP")
endif()

if (${WFMAXLINES})
    set_file_compile_options(waveFunction_m.f90 "-DWFMAXLINES=${WFMAXLINES}")
endif()

if (${WFMAXCOLS})
    set_file_compile_options(waveFunction_m.f90 "-DWFMAXCOLS=${WFMAXCOLS}")
endif()

if (${MAXDETSPERCSF})
    set_file_compile_options(multiDet_m.f90 "-DMAXDETSPERCSF=${MAXDETSPERCSF}")
endif()

target_link_libraries(src
        PUBLIC utils
        PUBLIC minimizers
        )

add_executable(amolqc amolqc_p.f90)

target_link_libraries(amolqc
        src)

add_fortran_library(inPsightsInterface STATIC
        amolqcInterface_m.f90
        )

target_link_libraries(inPsightsInterface
        src)

add_subdirectory(tests)
