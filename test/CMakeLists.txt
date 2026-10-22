# cmake driver file for testing

# Enable asserts
cet_enable_asserts()

# Add test items here

cet_test(ChannelOverlap_test
  SOURCE ChannelOverlap_test.cc
  LIBRARIES dunetrigger_Utilities
)