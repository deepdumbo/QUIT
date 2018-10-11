/*
 *  RFPulse.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 *  SPGR sequence with saturation pulse at different offsets
 *
 */

#include "RFPulse.h"

namespace QI {

RFPulse::RFPulse(const rapidjson::Value &json) {
  p1 = QI::GetMember(json, "p1").GetDouble();
  p2 = QI::GetMember(json, "p2").GetDouble();
}

rapidjson::Value RFPulse::toJSON(rapidjson::Document::AllocatorType &a) const {
  rapidjson::Value json_val(rapidjson::kObjectType);
  json_val.AddMember("p1", p1, a);
  json_val.AddMember("p2", p2, a);
  return json_val;
}

} // End namespace QI
