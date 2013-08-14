// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#include <Teuchos_TabularOutputter.hpp>

#include "MueLu_Level.hpp"

#include "MueLu_FactoryManagerBase.hpp"

namespace MueLu {

  RCP<Level> Level::Build() {
    RCP<Level> newLevel = rcp(new Level());

    // Copy 'keep' status of variables
    for (TwoKeyMap::const_iterator kt = map_.begin(); kt != map_.end(); kt++) {
      const FactoryBase* factory = kt->first;

      for (SubMap::const_iterator it = kt->second.begin(); it != kt->second.end(); it++) {
        const std::string& ename = it->first;

        if (IsKept(ename, factory, MueLu::Keep)) { // MueLu::Keep is the only flag propagated
          if (factory == NULL)                                  // TODO: Is this possible?? Throw exception. Not supposed to use the FactoryManager here.
            newLevel->Keep(ename, NoFactory::get());
          else
            newLevel->Keep(ename, factory);
        }
      }
    }

    return newLevel;
  }

  int Level::GetLevelID() const { return levelID_; }

  void Level::SetLevelID(int levelID) {
    if (levelID_ != -1 && levelID_ != levelID)
      GetOStream(Warnings1, 0) << "Warning: Level::SetLevelID(): Changing an already defined LevelID (previousID=" << levelID_ << ", newID=" << levelID << ")" << std::endl;

    levelID_ = levelID;
  }

  void Level::SetPreviousLevel(const RCP<Level> & previousLevel) {
    if (previousLevel_ != Teuchos::null && previousLevel_ != previousLevel)
      GetOStream(Warnings1, 0) << "Warning: Level::SetPreviousLevel(): PreviousLevel was already defined" << std::endl;

    previousLevel_ = previousLevel;
  }

  void Level::SetFactoryManager(const RCP<const FactoryManagerBase> & factoryManager) {
    factoryManager_ = factoryManager;
  }

  const RCP<const FactoryManagerBase> Level::GetFactoryManager() {
    return factoryManager_;
  }

  void Level::AddKeepFlag(const std::string& ename, const FactoryBase* factory, KeepType keep) {
    if (!IsKey(factory, ename)) {
      // If the entry does not exist, create it to store the keep flag
      Teuchos::RCP<MueLu::VariableContainer> newVar = Teuchos::rcp(new MueLu::VariableContainer);
      map_[factory][ename] = newVar;
    }
    // Set the flag
    map_[factory][ename]->AddKeepFlag(keep);
  }

  void Level::RemoveKeepFlag(const std::string& ename, const FactoryBase* factory, KeepType keep) {
    // No entry = nothing to do
    if (!IsKey(factory, ename))
      return;

    // Remove the flag
    Teuchos::RCP<MueLu::VariableContainer>& v = map_[factory][ename];
    v->RemoveKeepFlag(keep);

    // Remove data if no keep flag left and counter == 0
    if ((v->IsRequested() == false) && (v->GetKeepFlag() == 0)) {
      v = Teuchos::null; // free data

      map_[factory].erase(ename);
      if (map_.count(factory) == 0)
        map_.erase(factory);
    }
  }

  KeepType Level::GetKeepFlag(const std::string& ename, const FactoryBase* factory) const {
    if (!IsKey(factory,ename))
      return false;

    return Get(factory, ename)->GetKeepFlag();
  }

  void Level::Request(const FactoryBase& factory) {
    RequestMode prev = requestMode_;
    requestMode_ = REQUEST;
    factory.CallDeclareInput(*this);
    requestMode_ = prev;
  }

  void Level::Release(const FactoryBase& factory) {
    RequestMode prev = requestMode_;
    requestMode_ = RELEASE;
    factory.CallDeclareInput(*this);
    requestMode_ = prev;
  }

  void Level::DeclareInput(const std::string& ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
    if (requestMode_ == REQUEST) {
      try {
        Request(ename, factory, requestedBy);
      }
      catch(Exceptions::DependencyError &de) {
        std::string previousMsg(de.what());
        std::string msg = requestedBy->ShortClassName() + "::DeclareInput : (" + previousMsg + ") unable to find or generate requested data \""
                          + ename + "\"" + ((factory != NULL) ? " with generating factory " + factory->ShortClassName() + "." : ".");
        TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,msg);
        throw Exceptions::RuntimeError(msg);
      }
      catch(Exceptions::RuntimeError &rte) {
        std::string previousMsg(rte.what());
        std::string msg = previousMsg + "\n    during request for data \"" + ename + "\" by factory " + requestedBy->ShortClassName();
        TEUCHOS_TEST_FOR_EXCEPTION(true,Exceptions::RuntimeError,msg);
        throw Exceptions::RuntimeError(msg);
      }
    }
    else if (requestMode_ == RELEASE) {
      Release(ename, factory, requestedBy);
    }
    else
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareInput(): requestMode_ undefined.");
  }

  void Level::DeclareDependencies(const FactoryBase* factory, bool bRequestOnly, bool bReleaseOnly) { //TODO: replace bReleaseOnly, bReleaseOnly by one RequestMode enum
    if (bRequestOnly && bReleaseOnly)
      TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareDependencies(): Both bRequestOnly and bReleaseOnly set to true makes no sense.");

    if (requestMode_ == REQUEST) {

      if (bReleaseOnly == false) Request(*factory);

    } else if (requestMode_ == RELEASE) {

      if (bRequestOnly == false) Release(*factory);

    } else TEUCHOS_TEST_FOR_EXCEPTION(true, Exceptions::RuntimeError, "MueLu::Level::DeclareDependencies(): requestMode_ undefined.");
  }

  void Level::Request(const std::string& ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
    const FactoryBase* fac = GetFactory(ename, factory);

    // Historically, we called request for factory only if the factory has not been requested before and no data has
    // been generated by fact (independent of ename)
    //     bool test = (IsAvailableFactory(fac) == false && IsRequestedFactory(fac) == false);
    // However, this does not work for reuse of data, as some factories (like RAPFactory and EminPFactory) typically
    // have some generated data from previous setup (stencils, or initial prolonagor guess). Therefore, we remove
    // one check.
    // NOTE: It is possible that one might want to return this check when running in standalone mode, however
    // that would require determining when such situation happens which may be tricky.
    bool test = (IsRequestedFactory(fac) == false);

    // This request must be done before calling Request(*fac) to avoid circular dependency problems.
    if (!IsKey(fac, ename)) {
      Teuchos::RCP<MueLu::VariableContainer> newVar = Teuchos::rcp(new MueLu::VariableContainer);
      map_[fac][ename] = newVar;
    }

    Teuchos::RCP<MueLu::VariableContainer>& v = map_[fac][ename];
    v->Request(requestedBy);

    // The value of IsRequestedFactory(fac) is true, due to the above request.
    // That is why a temporary boolean "test" is used!
    TEUCHOS_TEST_FOR_EXCEPTION(IsRequestedFactory(fac) != true, Exceptions::RuntimeError, "Level::Request(ename, factory): internal logic error.");

    // Call Request for factory dependencies
    if (test)
      Request(*fac);
  }

  void Level::Release(const std::string& ename, const FactoryBase* factory, const FactoryBase* requestedBy) {
    const FactoryBase* fac = GetFactory(ename, factory);

    // Only a factory which has requested (fac,ename) is allowed to release it again.
    // Do not release data if it has not been requested by the factory "requestedBy"
    // Note: when data is released (fac,ename) depends on it often happened that some
    //       of this data has (recursively) been released too often
    if (IsRequestedBy(fac, ename, requestedBy)) {
      if (CountRequestedFactory(fac) == 1      &&     // check if factory fac is not requested by another factory
          IsAvailableFactory(fac)    == false) {      // check if Build function of factory fac has been called
        // In general all data (fac,ename) depends on is released with the Get calls in the Build
        // function of the generating factory fac.
        // Here we have to release the dependencies of some data that has been requested (by factory "requestedBy")
        // but the corresponding Build function of factory "fac" has never been called. Therefore the dependencies
        // have never been released. Do it now.
        Release(*fac);
      }

      TEUCHOS_TEST_FOR_EXCEPTION(!IsKey(fac,ename), Exceptions::RuntimeError, "\"" + ename + "\" not found. Do a request first.");

      Teuchos::RCP<MueLu::VariableContainer>& v = map_[fac][ename];
      v->Release(requestedBy);

      // Remove data if no keep flag left and counter == 0
      if ((v->IsRequested() == false) && (v->GetKeepFlag() == 0)) {
        v = Teuchos::null; // free data

        map_[fac].erase(ename);
        if (map_.count(fac) == 0)
          map_.erase(fac);
      }
    }
  }

  void Level::Clear() {
    // TODO: needs some love, ugly as it is
    // The ugliness is the fact that we restart both loops when we remove a single element
    bool wasRemoved;
    do {
      wasRemoved = false;
      for (TwoKeyMap::const_iterator kt = map_.begin(); kt != map_.end(); kt++) {
        const FactoryBase* factory = kt->first;

        for (SubMap::const_iterator it = kt->second.begin(); it != kt->second.end(); it++) {
          const std::string& ename = it->first;

          // We clear all the data that
          //   a) has not been requested
          //   b) is not being kept using NextRun (e.g., we clear out Final data)
          if (!IsKept(ename, factory, MueLu::NextRun)) {
            RemoveKeepFlag(ename, factory, MueLu::All); // will delete the data if counter == 0

            wasRemoved = true;
            break;
          }
        }

        if (wasRemoved)
          break;
      }

    } while (wasRemoved == true);
  }

  std::string Level::description() const {
    std::ostringstream out;
    out << BaseClass::description();
    out << "{ levelID = " << levelID_ << "}";
    return out.str();
  }

  void Level::print(Teuchos::FancyOStream& out, const VerbLevel verbLevel) const {
    RCP<Teuchos::FancyOStream> out0 = Teuchos::rcpFromRef(out);
    int previousSetting = out0->getOutputToRootOnly();
    out0->setShowProcRank(true);

    std::ostringstream ss;
    print(ss, verbLevel);

    out0->setOutputToRootOnly(-1);
    *out0 << ss.str();
    out0->setOutputToRootOnly(previousSetting);
    out0->setShowProcRank(false);
  }

  void Level::print(std::ostream& out, const VerbLevel verbLevel) const {
    out << "LevelID = " << GetLevelID() << std::endl;

    typedef Teuchos::TabularOutputter TTO;
    TTO outputter(out);
    outputter.pushFieldSpec("data name",                TTO::STRING, TTO::LEFT, TTO::GENERAL, 20);
    outputter.pushFieldSpec("gen. factory addr.",       TTO::STRING, TTO::LEFT, TTO::GENERAL, 18);
    outputter.pushFieldSpec("req",                      TTO::INT,    TTO::LEFT, TTO::GENERAL, 3);
    outputter.pushFieldSpec("keep",                     TTO::STRING, TTO::LEFT, TTO::GENERAL, 5);
    outputter.pushFieldSpec("type",                     TTO::STRING, TTO::LEFT, TTO::GENERAL, 15);
    outputter.pushFieldSpec("data",                     TTO::STRING, TTO::LEFT, TTO::GENERAL, 14);
    outputter.pushFieldSpec("req'd by",                 TTO::STRING, TTO::LEFT, TTO::GENERAL, 20);
    outputter.outputHeader();

    for (TwoKeyMap::const_iterator kt = map_.begin(); kt != map_.end(); kt++) {
      const FactoryBase* factory = kt->first;

      for (SubMap::const_iterator it = kt->second.begin(); it != kt->second.end(); it++) {
        const std::string& ename = it->first;

        outputter.outputField(ename);   // variable name

        // NOTE: We cannot dereference the factory pointer and call factory->description() as we do not know
        // if the factory still exist (the factory pointer is a raw pointer by design). Instead, the level
        // should store the factory description internally as a string for debugging purpose (and in debug mode only).
        //         // factory name
        //         std::stringstream ss1;
        //         ss1 << (*kt)->description();
        //         outputter.outputField((ss1.str()).substr(0,30));

        // factory ptr
        outputter.outputField(factory);

        int reqcount = NumRequests(factory, ename); // request counter
        outputter.outputField(reqcount);

        KeepType keepType = GetKeepFlag(ename, factory);
        if (keepType != 0) {
          std::stringstream ss;
          if (keepType & MueLu::UserData) { ss << "User";  }
          if (keepType & MueLu::Keep)     { ss << "Keep";  }
          if (keepType & MueLu::Final)    { ss << "Final"; }
          outputter.outputField(ss.str());
        } else {
          outputter.outputField("No");
        }

        if (IsAvailable(ename, factory)) {
          std::string strType = it->second->GetData().getAny(true).typeName();

          if (strType.find("Xpetra::Matrix") != std::string::npos) {
            outputter.outputField("Matrix" );
            outputter.outputField("available");
          } else if (strType.find("Xpetra::MultiVector") != std::string::npos) {
            outputter.outputField("Vector");
            outputter.outputField("available");
          } else if (strType.find("Xpetra::Map") != std::string::npos) {
            outputter.outputField("Map");
            outputter.outputField("available");
          } else if (strType.find("MueLu::SmootherBase") != std::string::npos) {
            outputter.outputField("SmootherBase");
            outputter.outputField("available");
          } else if (strType.find("MueLu::Aggregates") != std::string::npos) {
            outputter.outputField("Aggregates");
            outputter.outputField("available");
          } else if (strType.find("MueLu::AmalgamationInfo") != std::string::npos) {
            outputter.outputField("AmalgamationInfo");
            outputter.outputField("available");
          } else if (strType.find("MueLu::Graph") != std::string::npos) {
            outputter.outputField("Graph");
            outputter.outputField("available");
          } else if (strType.find("MueLu::Constraint") != std::string::npos) {
            outputter.outputField("Constraint");
            outputter.outputField("available");
          } else if (strType == "int") {
            outputter.outputField(strType);
            outputter.outputField(Teuchos::getValue<int>(it->second->GetData()));
          } else if (strType == "double") {
            outputter.outputField(strType);
            outputter.outputField(Teuchos::getValue<double>(it->second->GetData()));
          } else if (strType == "string") {
            outputter.outputField(strType);
            outputter.outputField(Teuchos::getValue<std::string>(it->second->GetData()));
          } else {
            outputter.outputField(strType);
            outputter.outputField("available");
          }
        } else {
          outputter.outputField("unknown");
          outputter.outputField("not available");
        }

        typedef VariableContainer::request_container container_type;
        const container_type& requestedBy = it->second->Requests();
        std::ostringstream ss;
        for (container_type::const_iterator ct = requestedBy.begin(); ct != requestedBy.end(); ct++) {
          if (ct != requestedBy.begin())    ss << ",";
                                            ss << ct->first;
          if (ct->second > 1)               ss << "(" << ct->second << ")";
        }
        outputter.outputField(ss.str());

        outputter.nextRow();
      }
    } // for (TwoKeyMap::const_iterator kt = map_.begin(); kt != map_.end(); kt++) {
  }

#if defined(HAVE_MUELU_BOOST) && defined(BOOST_VERSION) && (BOOST_VERSION >= 104400)
    void Level::UpdateGraph(std::map<const FactoryBase*, BoostVertex>&                   vindices,
                            std::map<std::pair<BoostVertex, BoostVertex>, std::string>&  edges,
                            BoostProperties&                                             dp,
                            BoostGraph&                                                  graph) const {
      size_t vind = vindices.size();

      for (TwoKeyMap::const_iterator it1 = map_.begin(); it1 != map_.end(); it1++) {
        if (vindices.find(it1->first) == vindices.end()) {
          BoostVertex boost_vertex = boost::add_vertex(graph);
          boost::put("label", dp, boost_vertex, it1->first->description());
          vindices[it1->first] = vind++;
        }

        for (SubMap::const_iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++) {
          const VariableContainer::request_container& requests = it2->second->Requests();
          for (VariableContainer::request_container::const_iterator rit = requests.begin(); rit != requests.end(); rit++) {
            if (vindices.find(rit->first) == vindices.end()) {
              // requested by factory which is unknown
              BoostVertex boost_vertex = boost::add_vertex(graph);
              boost::put("label", dp, boost_vertex, rit->first->description());
              vindices[rit->first] = vind++;
            }

            edges[std::pair<BoostVertex,BoostVertex>(vindices[rit->first], vindices[it1->first])] =  it2->first;
          }
        }
      }
    }
#endif

  // JG Note: should the option IgnoreUserData() moved to the Factory interface or on the specific factories that are using this option? It would simplify the level class.
  const FactoryBase* Level::GetFactory(const std::string& ename, const FactoryBase* factory) const {
    if (factory != NULL)
      return factory;

    // If IgnoreUserData == false and if variable "ename" generated by NoFactory is provided by the user (MueLu::UserData),
    // use user-provided data by default without querying the FactoryManager.
    // When FactoryManager == null, we consider that IgnoreUserData == false.
    if ((factoryManager_ == Teuchos::null || factoryManager_->IgnoreUserData() == false) &&
        (IsAvailable(ename, NoFactory::get()) && IsKept(ename, NoFactory::get(), MueLu::UserData))) {
      return NoFactory::get();
    }

    // Query factory manager
    TEUCHOS_TEST_FOR_EXCEPTION(factoryManager_ == null, Exceptions::RuntimeError, "MueLu::Level("<< levelID_ << ")::GetFactory(" << ename << ", " << factory << "): No FactoryManager");
    const FactoryBase* fac = factoryManager_->GetFactory(ename).get();
    TEUCHOS_TEST_FOR_EXCEPTION(fac == NULL, Exceptions::RuntimeError, "MueLu::Level("<< levelID_ << ")::GetFactory(" << ename << ", " << factory << "): Default factory returned by FactoryManager cannot be NULL");
    return fac;
  }

  Level::RequestMode Level::requestMode_ = UNDEF;

} //namespace MueLu

//TODO: Caps should not matter
