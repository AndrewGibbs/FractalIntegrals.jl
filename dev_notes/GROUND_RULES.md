# Ground rules for AI agents

These come from Andrew and apply to any agent working on this repo. They replace the shorter list in `next_steps.md`. They exist because **the humans must keep understanding the codebase**: AI should not outrun them with ideas, derivations or code.

## What may change

1. **New features live in new subfolders of `src/`.** Keep AI contributions separate from existing code.
2. **Change existing code only when one of these applies, and keep the change minimal:**
   - a new feature needs it for compatibility;
   - it fixes a bug. Show Andrew the bug and get his approval first, then fix it with a test, one small fix at a time;
   - it improves performance. Work on one function at a time, keep the change only if a benchmark shows a real gain, and show Andrew the diff. Do not restructure his code for small gains, because he must still be able to read and understand it.
3. **Record every change to existing code** in `dev_notes/changes_to_existing_code.md`.
4. **Ask before adding any dependency.** Import time is already a concern.
5. **Do not edit docs or docstrings for existing code** unless Andrew asks. He and Josh have not yet decided whether AI should write them.

## What every change needs

6. **Unit tests in `test/`.** Testing is fully delegated to agents, who decide what is appropriate. All tests must pass before a change is finalised.
7. **New features need thorough documentation that a human reads.** Store the mathematical derivation in `dev_notes/derivations/`.

## Who decides

8. **Only humans authorise official version releases.**
9. **Commit or push only when asked.**

## How to work

- Write the smallest thing that works, and match the surrounding style (Unicode maths names, comment density, naming).
- Work at a pace the humans can review. Explain derivations and design choices before building on them, and do not push on to see the end result first.
- Do not start planned features until asked.
- Treat a finding as unconfirmed until it has been run, and report how much it matters as well as whether it is true.
- Andrew and Josh work on separate branches and merge into `main` when features are complete. Check that work does not overlap with the other person's. Josh leads the FFT-Toeplitz work on his branch.
