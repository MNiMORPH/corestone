# <decision being made>

Copy this to `design/<short-name>.md` and fill it in **before** writing the
implementation.

## The question

What has to be decided, in one or two sentences. If it is not actually a
decision -- if there is one obvious way -- do not write a design doc; just
write the code.

## Options, each with its full cost

For each: what it does, and its whole blast radius -- the tests it forces, the
docs it ripples into, the measured results it invalidates. A choice should
carry its true cost at the moment it is made.

## The probe

The smallest runnable thing that settles this, and what it would show either
way. Write it in `prototypes/`, run it, and paste the output here. A design
doc without a probe is a guess.

## Decision, and why

Including which option was rejected and what would change the answer.

## Parameters chosen here

Every threshold, cut-off, filter, or default this decision introduces, with its
source. Anything not asked for is a proposal, flagged as such.
