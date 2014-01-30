Module Timing

Type, Public :: Timer
	Real*8 :: delta, elapsed
    Real*4 :: t1

	Contains
	Procedure :: Init  => Initialize_Timer
	Procedure :: startclock
	Procedure :: stopclock
	Procedure :: increment
End Type Timer

Contains

Subroutine Initialize_Timer(self)
	Implicit None
	Class(Timer) :: self
	self%elapsed = 0.0d0
	self%t1 = 0.0d0
	self%delta = 0.0d0
End Subroutine Initialize_Timer

Subroutine Startclock(self)
	Implicit None
	Class(Timer) :: self
	self%t1 = secnds(0.0)
End Subroutine Startclock

Subroutine Stopclock(self)
	Implicit None
	Class(Timer) :: self
	self%delta = secnds(self%t1)
End Subroutine Stopclock

Subroutine Increment(self)
	Implicit None
	Class(Timer) :: self
	Call self%stopclock()
	self%elapsed = self%elapsed+self%delta
End Subroutine Increment
End Module Timing
