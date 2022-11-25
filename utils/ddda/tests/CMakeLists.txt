add_pFUnit_ctest(ddda_test
        TEST_SOURCES boltzmann_tm.pf decorrelation_tm.pf statistic_tm.pf
        LINK_LIBRARIES basics ddda optimizers generic
        LABELS utils moduletest Amolqc pFUnit)
