/*
 *  MultiEchoSequence.cpp
 *
 *  Copyright (c) 2018 Tobias Wood.
 *
 *  This Source Code Form is subject to the terms of the Mozilla Public
 *  License, v. 2.0. If a copy of the MPL was not distributed with this
 *  file, You can obtain one at http://mozilla.org/MPL/2.0/.
 *
 */

#include "MultiEchoSequence.h"
#include "Macro.h"

namespace QI {

/*
 * Base
 */
Eigen::Index MultiEchoBase::size() const {
    return TE.rows();
}

/*
 * Regularly spaced sequence
 */
MultiEchoSequence::MultiEchoSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = GetMember(json, "TR").GetDouble();
    TE1 = GetMember(json, "TE1").GetDouble();
    ESP = GetMember(json, "ESP").GetDouble();
    ETL = GetMember(json, "ETL").GetInt();
    TE = Eigen::ArrayXd::LinSpaced(ETL, TE1, TE1+ESP*(ETL - 1));
}

rapidjson::Value MultiEchoSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("TE1", TE1, a);
    json.AddMember("ESP", ESP, a);
    json.AddMember("ETL", ETL, a);
    return json;
}

/*
 * Irregularly spaced sequence
 */
MultiEchoFlexSequence::MultiEchoFlexSequence(const rapidjson::Value &json) {
    if (json.IsNull()) QI_FAIL("Could not read sequence: " << name());
    TR = GetMember(json, "TR").GetDouble();
    TE = ArrayFromJSON(json, "TE");
}

rapidjson::Value MultiEchoFlexSequence::toJSON(rapidjson::Document::AllocatorType &a) const {
    rapidjson::Value json(rapidjson::kObjectType);
    json.AddMember("TR", TR, a);
    json.AddMember("TE", ArrayToJSON(TE, a), a);
    return json;
}

} // End namespace QI
