#ifndef DUNETRIGGER_FDTPC_TPCTRIGGERSIM_VERBOSITY_HH
#define DUNETRIGGER_FDTPC_TPCTRIGGERSIM_VERBOSITY_HH

namespace dunetrigger::fdtpc {
enum Verbosity {
  // these would all implicitly have these values
  // but it's best to explicitly define things
  // for anyone in the future looking
  kQuiet = 0,
  kInfo = 1,
  kDebug = 2,
  kVerbose = 3
};
} // namespace dunetrigger::fdtpc

#endif
