# C++ Coding Style and Conventions for FDPS

## General Principles
- Follow modern C++ best practices
- Respect clang-tidy profile and treat warnings as errors
- Use Serena MCP for understanding codebase structure and reuse
- Prefer simplest correct design addressing root causes
- Write modular code with single responsibility
- Practice test-driven development (Given, When, Then style)

## Naming Conventions
- Types and enums: PascalCase
- Functions and methods: camelCase
- Variables and data members: snake_case
- Constants: SCREAMING_SNAKE_CASE or kCamelCase
- Filenames: snake_case.hpp or snake_case.cpp
- Accurate filenames reflecting primary purpose

## Code Structure
- No hard-coded strings: use constexpr constants, enum class, string_view
- No magic numbers: use named constexpr constants or enum class
- Floating point: no direct equality comparison, use tolerant comparisons, prefer double, guard with std::isfinite
- Comments: explain reasoning and non-obvious parts, not restating code
- Documentation: for public APIs, update docs, include change log and root-cause for fixes

## Modern C++ Features
- Prefer RAII, std::unique_ptr, std::shared_ptr, std::optional, std::variant
- Ranges and algorithms libraries
- Avoid raw new/delete, implicit conversions
- Small, pure functions, isolate I/O/side effects
- enum class over unscoped enums
- constexpr over preprocessor defines

## Macros (Limited Use)
- Only when no safer alternative
- SCREAMING_SNAKE_CASE with project prefix
- No shadowing standard/third-party identifiers
- No leading/double underscores
- Parenthesize parameters and expressions
- No side effects, no multiple evaluation
- For multi-statement: do { ... } while (0)
- Document purpose, inputs, example, why macro needed
- Include guards: unique project-scoped name
- Audit periodically for issues

## Architecture
- No new compatibility layers
- Remove legacy code in same change when practical
- Confirm file placement with Serena MCP
- Avoid repetition: extract common logic

## Testing
- Unit tests first for new features/bug fixes
- Cover boundary conditions, error cases, invariants
- GoogleTest framework implied

## Repository Hygiene
- Remove .bak, .orig, swap files, temp logs, test screenshots
- Accurate .gitignore

## Change Completion Checklist
- Tests green, written first
- No hard-coded strings/magic numbers
- Builds cleanly with linting/warnings as errors
- Folder/naming conform to rules (verify with Serena MCP)
- Docs updated, root-cause noted for fixes
- Legacy code removed, no new compatibility layers
- CI passes on supported toolchains