#ifndef DUNETRIGGER_TRIGGERSIM_GETMANYBYREGEXTAG_HH
#define DUNETRIGGER_TRIGGERSIM_GETMANYBYREGEXTAG_HH

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Selector.h"
#include "canvas/Persistency/Provenance/BranchDescription.h"
#include "canvas/Utilities/InputTag.h"

#include <regex>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <vector>

namespace dunetrigger {

  template<typename T>
  std::vector<art::Handle<T>>
  getManyByRegexTag(art::Event& e, const art::InputTag& tag)
  {
    std::regex instance_regex(!tag.instance().empty() ? tag.instance() : ".*");
    std::regex label_regex   (!tag.label()   .empty() ? tag.label()    : ".*");
    std::regex process_regex (!tag.process() .empty() ? tag.process()  : ".*");

    art::SelectorByFunction selector(
      [instance_regex, label_regex, process_regex](art::BranchDescription const& p) {
        return std::regex_match(p.inputTag().label(),    label_regex)
             & std::regex_match(p.inputTag().instance(), instance_regex)
             & std::regex_match(p.inputTag().process(),  process_regex);
      },
      "InputTag Regex Selector"
    );

    auto handles = e.getMany<T>(selector);
    if (handles.empty()) {
      throw std::runtime_error(
        "No " + std::string(typeid(T).name()) +
        " collections matching " + tag.instance() + "_" + tag.label());
    }
    return handles;
  }

} // namespace dunetrigger

#endif // DUNETRIGGER_TRIGGERSIM_GETMANYBYREGEXTAG_HH
