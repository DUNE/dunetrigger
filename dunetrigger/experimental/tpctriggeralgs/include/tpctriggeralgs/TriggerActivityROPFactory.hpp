/* @file: TriggerActivityROPFactory.hpp
 *
 * This is part of the DUNE DAQ Application Framework, copyright 2023.
 * Licensing/copyright details are in the COPYING file that you should have
 * received with this code.
 */

#ifndef tpctriggeralgs_TRIGGER_ACTIVITY_ROP_FACTORY_HPP_
#define tpctriggeralgs_TRIGGER_ACTIVITY_ROP_FACTORY_HPP_

#include "dunetrigger/experimental/tpctriggeralgs/include/tpctriggeralgs/TriggerActivityROPMaker.hpp"
#include "AbstractFactory.hpp"

#define REGISTER_TRIGGER_ACTIVITY_ROP_MAKER(tam_name, tam_class)                                                                                      \
  static struct tam_class##Registrar {                                                                                                            \
    tam_class##Registrar() {                                                                                                                      \
      TriggerActivityROPFactory::register_creator(tam_name, []() -> std::unique_ptr<TriggerActivityROPMaker> {return std::make_unique<tam_class>();});   \
    }                                                                                                                                             \
  } tam_class##_registrar;

namespace tpctriggeralgs {

class TriggerActivityROPFactory : public AbstractFactory<TriggerActivityROPMaker> {};

} /* namespace tpctriggeralgs */

#endif // tpctriggeralgs_TRIGGER_ACTIVITY_ROP_FACTORY_HPP_
