#ifndef MUELU_MUTUALLYEXCLUSIVETIME_HPP
#define MUELU_MUTUALLYEXCLUSIVETIME_HPP

#include <string>
#include <stack>
#include <Teuchos_Time.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "MueLu_Exceptions.hpp"
#include "MueLu_BaseClass.hpp"

namespace MueLu {

  //! This class wraps a Teuchos::Time and maintains a mutually exclusive property between wrapped timers.
  //! When a MutuallyExclusiveTime is running, other timers are not running.
  //! Timers have three state: running, stopped or paused to enforce the mutually exclusive property.
  //! When the running timer is stopped, the last active timer is restarted. A stack of timers is used internally to support this functionality. 
  //! This class is useful to exclude from a timer the execution time of a subroutine.
  //!
  //! Example:
  //!
  //! Note: Only one timer can be active at a time but all timers can be inactive at the same time. Timers cannot be destroyed when they are in 'paused'.

  //TODO: inheritence from PerformanceMonitorBase<Time> ?

  template<class TagName> //! The template parameter of this class can be used to define several set of mutually exclusive timer.
  class MutuallyExclusiveTime : public BaseClass {
    
  public:
    
    MutuallyExclusiveTime(const std::string &name, bool start=false)
      : timer_(rcp(new Teuchos::Time(name, false))),  // second argument is false in any case, because if start==true, 
                                                      // timer has to be started by MutuallyExclusiveTime::start() instead of Teuchos::Time::start().
	isPaused_(false)
    { 
      if (start == true) timer_->start();
    }
    
    ~MutuallyExclusiveTime() {
      // This timer can only be destroyed if it is not in the stack
      if (isPaused()) {
	// error message because cannot throw an exception in destructor
	GetOStream(Errors, 0) << "MutuallyExclusiveTime::~MutuallyExclusiveTime(): Error: destructor called on a paused timer." << std::endl;
	//TODO: Even if timing results will be wrong, the timer can be removed from the stack to avoid a segmentation fault.
      }
      
      stop(); // if isRunning(), remove from the stack, resume previous timer
    }
    
    //! Starts the timer. If a MutuallyExclusiveTime is running, it will be stopped.
    //! Precondition: timer is not already paused
    //! Postcondition: timer is running. Other MutuallyExclusiveTime are paused or stop.
    void start(bool reset=false) {
      TEUCHOS_TEST_FOR_EXCEPTION(isPaused(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::start(): timer is paused. Use resume().");
      
      if (isRunning()) { return; } // If timer is already running, do not pause/push-in-the-stack/start the timer. 
                                   // Otherwise, something bad will happen when this.stop() will be called
      
      // pause currently running timer
      if (!timerStack_.empty()) {
	timerStack_.top()->pause();
      }
      
      // start this timer
      timer_->start(reset);
      timerStack_.push(this);
    }
    
    // {@ Functions that can only be called on the most recent timer (= running timer or last paused timer)
    
    //!	Stops the timer. The previous MutuallyExclusiveTime that has been paused when this timer was started will be resumed.
    // stop() can be called on an already stopped timer or on the currently running timer
    double stop() {
      TEUCHOS_TEST_FOR_EXCEPTION(isPaused(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::start(): timer is paused. Use resume().");
      if (!isRunning()) { return timer_->stop(); } // stop() can be called on stopped timer
      
      // Here, timer is running, so it is the head of the stack
      TopOfTheStack();
      
      timerStack_.pop();
      double r = timer_->stop();

      if (!timerStack_.empty()) {
	timerStack_.top()->resume();
      }
      
      return r;
    }
    
    //! Pause running timer. Used internally by start().
    void pause() {
      if (isPaused()) // calling twice pause() is allowed
        return;

      TopOfTheStack();
      
      timer_->stop();
      isPaused_ = true;
    }
    
    //! Resume paused timer. Used internally by stop()
    //! Precondition: timer is at the top of the stack
    //! Timer is not reset
    void resume() {
      TopOfTheStack();
      
      // no 'shortcut' test necessary: 
      // - if timer is stop, it is in pause (cannot be stop and not in pause because this timer is the head of the stack).
      // - if timer is running, nothing is changed by this function.

      timer_->start(false);
      isPaused_ = false;
    }

    // @}    

    //@{

    bool isRunning() {
      if (timer_->isRunning()) {
	TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.top() != this, Exceptions::RuntimeError, 
				   "MueLu::MutuallyExclusiveTime::isRunning(): this timer is active so it is supposed to be the head of the stack");
      }
      return timer_->isRunning();
    }
    
    bool isPaused() {
      TEUCHOS_TEST_FOR_EXCEPTION(isPaused_ && timer_->isRunning(), Exceptions::RuntimeError, "");
      return isPaused_;
    }

    //@}

    //! Return a new MutuallyExclusiveTime that is register with the Teuchos::TimeMonitor (for timer summary)
    // Note: this function is provided by the timer class, not by a monitor (!= Teuchos)
    static RCP<MutuallyExclusiveTime<TagName> > getNewTimer(const std::string& name) {
      RCP<MutuallyExclusiveTime<TagName> > timer = rcp(new MutuallyExclusiveTime<TagName>(Teuchos::TimeMonitor::getNewTimer(name)));
      return timer;
    }

    //! Increment the number of times this timer has been called. 
    void incrementNumCalls() { timer_->incrementNumCalls(); }

  private:
    
    // This constructor is not public because I'm concerned that users will used Teuchos::Time::start()/stop() 
    // instead of MutuallyExclusiveTime::start()/stop() if they have access to the underlying Teuchos::Time object.
    MutuallyExclusiveTime(RCP<Teuchos::Time> timer)
      : timer_(timer), isPaused_(false)
    { }
    
    // MutuallyExclusiveTime() { }
    
    RCP<Teuchos::Time> timer_; // using an RCP allows to use Teuchos::TimeMonitor to keep track of the timer.
    bool isPaused_;

    //! Stack of started timers (active or paused timers).
    // - empty when no active timer
    // - head is the active timer
    // - other timers are timers paused to enforce the mutually exclusive property of the timer set.
    static std::stack<MutuallyExclusiveTime<TagName>*> timerStack_;

    //! Check if 'this' is the head of the stack.
    void TopOfTheStack() {
      TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.empty(), Exceptions::RuntimeError, "MueLu::MutuallyExclusiveTime::TopOfTheStack(): timer is not the head of the stack (stack is empty).");
      TEUCHOS_TEST_FOR_EXCEPTION(timerStack_.top() != this, Exceptions::RuntimeError,  "MueLu::MutuallyExclusiveTime::TopOfTheStack(): timer is not the head of the stack.");
      TEUCHOS_TEST_FOR_EXCEPTION(!(isRunning() || isPaused()), Exceptions::RuntimeError,  "MueLu::MutuallyExclusiveTime::TopOfTheStack(): head of the stack timer is neither active nor paused.");
    }

  //TODO: test integrity of the stack:
  // Head = running or paused
  // Other timers of the stack = paused
  
  };

} // namespace MueLu

#endif // MUELU_MUTUALLYEXCLUSIVETIME_HPP


