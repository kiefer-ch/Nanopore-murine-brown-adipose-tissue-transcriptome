rule packrat_init:
    run:
        shell("""
            R -e 'source("packrat/init.R")'
        """)
        shell("R -e 'packrat::restore()'")
