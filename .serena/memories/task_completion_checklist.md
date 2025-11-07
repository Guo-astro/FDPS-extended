# Task Completion Checklist for FDPS

Before marking a task as complete, ensure:

## Code Quality
- [ ] Mapped structure with Serena MCP, reused existing facilities
- [ ] Removed compatibility layers and legacy code in this change
- [ ] Eliminated hard-coded strings and magic numbers
- [ ] Conformed to naming and folder conventions (verified with Serena MCP)
- [ ] Handled floating point with tolerant comparisons and std::isfinite
- [ ] Written comments explaining reasoning, documented root cause for fixes
- [ ] Followed macro rules (limited use, proper naming, documentation)

## Testing
- [ ] Written tests first in Given, When, Then style
- [ ] Tests pass for all new/modified behavior
- [ ] Covered boundary conditions, error cases, invariants

## Build and Linting
- [ ] Builds cleanly with clang-tidy warnings as errors
- [ ] Compiled with -Werror -Wall -Wextra -Wpedantic
- [ ] Passed sanitizers (address, undefined behavior) in CI if feasible

## Repository
- [ ] Cleaned workspace of artifacts (.bak, .orig, .tmp, .log files)
- [ ] Accurate .gitignore maintained
- [ ] No misplaced files (validated against allowed folders)

## CI/CD
- [ ] Passes continuous integration across supported toolchains
- [ ] GoogleTest output CI-friendly

## Documentation
- [ ] Updated documentation for changes
- [ ] Included concise change log
- [ ] Root-cause summary and best practice explanation for fixes

## Commands to Run
- Build: `cmake --build .` or `ninja`
- Test: `ctest` (if configured)
- Lint: `clang-tidy <files>` with project profile
- Format: (if clang-format configured)
- Clean: Remove build artifacts before commit