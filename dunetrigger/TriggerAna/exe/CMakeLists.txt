cet_test(test_tuple_writing
  SOURCE
    test_tuple_writing.cxx
  LIBRARIES
    larcore::headers
    Boost::headers
    ROOT::Tree
  NO_AUTO
)

install_source()
