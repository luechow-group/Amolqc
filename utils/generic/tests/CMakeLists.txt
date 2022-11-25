add_pFUnit_ctest(generic_test
        TEST_SOURCES blockAllocator_tm.pf genericFilter_tm.pf plateau_tm.pf
        LINK_LIBRARIES basics generic
        LABELS utils moduletest Amolqc pFUnit)
