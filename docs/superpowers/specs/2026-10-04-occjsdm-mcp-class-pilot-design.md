# occJSDM MCP class pilot design

Status: draft for Doug's review, 4 October 2026.

## Intended outcome

For the November class in Indonesia, students can follow the existing occJSDM lessons and use an AI assistant to run a supported occJSDM analysis on their own data. The preferred delivery is a Posit Cloud RStudio project with Posit Assistant and each student's OpenRouter account, so students do not need to install R or maintain a paid model subscription. This design extends the existing lessons rather than replacing them.

## What is known and what remains an assumption

Doug confirmed the class is in Indonesia, free OpenRouter accounts are acceptable, and the project should use Posit Cloud and Posit Assistant if feasible. The current Posit Assistant documentation describes OpenRouter as a supported provider and allows MCP server configuration. Posit Cloud documentation says Posit Assistant requires a Posit AI subscription. It is not yet confirmed whether a student can use an individual OpenRouter key inside Posit Assistant running in Posit Cloud, or whether the Cloud environment permits the required MCP process and dependencies. These are feasibility questions for the first pilot, not settled design facts.

The MCP server would run in each student's Posit Cloud project container. In this sense it is local to the RStudio session, even though the project itself runs remotely in Posit Cloud. Posit Assistant supports project-level local MCP server configuration, but project settings apply only in a trusted workspace. The class setup must show students how to trust the prepared project and must check whether Posit Cloud allows the needed local subprocess. The MCP server is an R-facing interface to selected occJSDM operations, not a new statistical implementation and not arbitrary R code execution.

## Proposed architecture

Use the existing occJSDM R package as the statistical engine. Add an optional MCP component in the same repository, using the R MCP implementation `mcptools` where its current APIs and Posit Assistant compatibility support the required workflow. Keep MCP dependencies out of the core package's required installation path. The MCP component should expose a small set of named, documented tools and call existing package functions.

Students open a prepared Posit Cloud project containing the optional MCP component, a small example dataset, and project-level MCP configuration. Each student configures their own OpenRouter key through Posit Assistant's provider setup and credential mechanism; no key is committed to the repository or project template. Posit Assistant calls the tools through MCP; the tools execute R code in that student's Cloud project and save fit artifacts there.

The assistant should receive concise validation reports and analysis summaries. The MCP tools should not transmit full datasets, arbitrary file contents, or full posterior arrays to the model by default. Model requests, tool descriptions, summaries, and any data sent to OpenRouter are subject to the provider's terms and the student's account settings; the class guide must make this clear before students use their own data.

## First class workflow

1. Check an occJSDM input object for required `info` and `OTU` structure, identify supported model structure, and report actionable issues before fitting.
2. Fit a non-spatial occJSDM model with explicit, class-safe MCMC settings and save the fit in the project.
3. Return convergence diagnostics and a compact, clearly labelled summary for discussion in Posit Assistant.

The first version targets an R object already shaped for `runOccJSDM()`, stored as an `.rds` file in the student's project. The pilot must verify this workflow with both the bundled example and a separately prepared user dataset. Additional raw-file formats, automatic covariate formula creation, maps, predictions, and broad tool catalogs are outside the first version unless the pilot shows they are essential to teach the workflow.

Do not expose spatial fitting in the initial MCP workflow. The package's spatial length-scale update is known to score the wrong density, so spatial fits cannot currently be presented as reliable. Do not expose arbitrary R evaluation or unrestricted filesystem access.

## Class and lesson integration

Preserve the existing lesson sequence and add a short, clearly optional AI-assisted analysis extension after the existing fitting and diagnostics lessons. Reuse existing teaching datasets and explanations where they fit. The extension should explain the MCP tool workflow, accepted input object, fit settings, how to read convergence output, and how to inspect the saved R result directly. It should show that the assistant can help operate and explain the package while the statistical estimates still come from occJSDM.

Provide a non-AI route through the existing lessons and R functions. It should remain possible to teach if students cannot access a free model, encounter rate limits, or cannot use the MCP integration.

## Feasibility gate and success criteria

The first deliverable is a disposable technical pilot, not package integration. It succeeds only if a student-like Posit Cloud account can enable Posit Assistant, configure OpenRouter using an individual key, connect to an MCP server running in the project, invoke a tool, and receive the expected result. Verify a complete small non-spatial fit and diagnostic summary in the Cloud project's resource and runtime limits. Record the versions, exact setup steps, model identifier, and any limitations. If any required connection is unavailable, revise the client or hosting architecture before integrating it into the package or class project.

The class-ready workflow is successful when a student can start from the shared project, provide an eligible RDS input, get validation feedback, run the supported fit, retrieve diagnostics, and locate the saved result without editing MCP configuration by hand. The same lesson must remain teachable through the non-AI route.

## Scope boundaries

- Optional MCP support lives alongside occJSDM and does not make MCP or `mcptools` mandatory for ordinary package users.
- The AI client is Posit Assistant in RStudio on Posit Cloud, contingent on the feasibility gate.
- The model provider is OpenRouter using student-owned API keys and currently available free models that support tool calls.
- The first workflow supports an occJSDM-ready RDS object and non-spatial fitting plus diagnostics.
- The MCP layer uses allowlisted tools with constrained inputs, bounded runtime settings, and project-scoped file access.
- Existing lessons remain the teaching foundation; add one optional lesson extension and a guided setup/checklist.
- No central hosting service, shared class API key, paid model subscription, arbitrary R execution, spatial fitting, or broad analysis assistant is part of this pilot.

## Risks and mitigations

Posit Cloud may require Posit AI Pass even when OpenRouter is configured. Resolve this in the feasibility gate before relying on it for the class. Project MCP settings may require students to trust the project, and the Cloud environment may restrict launching local subprocesses; verify both. Free OpenRouter models can change, be rate-limited, or fail to call tools consistently; pin a verified model for the lesson where possible and retain the non-AI route. Cloud CPU, memory, or runtime limits may make MCMC fits impractical; use an explicit small teaching configuration and validate runtime, and retain precomputed results for the lesson. Students may use data they should not send to a third-party model; keep raw data out of tool responses and explain what content the assistant receives. Input format mistakes can block novices; start from an RDS object in the package's native input shape and give a worked conversion example in the lesson.

## Decisions needed after the pilot

If the Posit Cloud and OpenRouter integration works, confirm the exact model, resource settings, supported input structure, data handling language, and class project distribution method before implementation. If it does not work, choose between RStudio Desktop plus a local MCP server, another MCP-capable client, or a non-MCP Posit Cloud workflow. No fallback is selected until the pilot identifies the actual blocker.

## References

- Posit Assistant provider setup, including OpenRouter: https://assistant.posit.co/docs/getting-started/providers/
- Posit Assistant MCP configuration, including local stdio servers and trusted project settings: https://assistant.posit.co/docs/reference/mcp-servers/
- Posit Assistant configuration reference, including project-level MCP settings: https://assistant.posit.co/docs/reference/config-file/
- Posit Cloud AI setup and its current Posit AI subscription note: https://docs.posit.co/cloud/get_started/index.html
- `mcptools` R server implementation and Posit Connect deployment route: https://github.com/posit-dev/mcptools/blob/main/R/server.R
- Posit Cloud project containers and resource limits: https://docs.posit.co/cloud/guide/projects/index.html
