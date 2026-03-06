cet_make_exec(soa_demo
  SOURCE
    soa_demo.cxx
  LIBRARIES
    larcore::headers
    Boost::headers
    ROOT::Tree
)

install_source()
