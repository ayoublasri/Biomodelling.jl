# Security policy

## Scope

Biomodelling.jl is a simulation library. It does not run a service, open a socket,
authenticate anyone, or handle credentials, so most classes of vulnerability do not
apply. What is in scope:

- the Julia package in `src/` and its dependencies;
- the GitHub Actions workflows in `.github/workflows/`;
- the Node tooling in `video/`, which builds a promotional video and is never
  executed by the library or by users of it.

`video/` is a build-time tool. Nothing it depends on is shipped to users or reachable
at run time, so an advisory against it is a housekeeping matter rather than a risk to
anyone running the package.

## Reporting

Report suspected vulnerabilities privately through GitHub's
[security advisory form](https://github.com/ayoublasri/Biomodelling.jl/security/advisories/new)
rather than in a public issue. Please include a reproduction and the versions
involved. Expect an acknowledgement within a week.

## Supported versions

Fixes land on the default branch and in the next release. Only the latest release
is supported.

## What the project does to stay clean

- Dependabot watches the npm tree and the Actions workflows weekly
  (`.github/dependabot.yml`).
- Every workflow declares least-privilege `permissions:`.
- CI fails if `npm audit` reports a moderate or worse advisory in `video/`.
- No secrets, keys or font binaries are committed; the build fetches what it needs.
