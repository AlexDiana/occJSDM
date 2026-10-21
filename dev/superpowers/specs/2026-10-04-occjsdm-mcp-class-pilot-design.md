# occJSDM MCP class pilot design

Status: Doug approved the implementation plan on 5 October 2026. First written 4 October 2026; updated 5 October 2026 with the Paper2Agent approach, workshop model fallback, lesson-guidance goal and model-cost comparison. Local baseline preparation and an unscored direct-R Cloud smoke test have begun; the small fit/save/diagnostics/reload sequence completed with setup assistance. Local MCP conversion and independent verification have passed on macOS arm64, including real stdio calls with a separately restored native runtime. Full Cloud feasibility testing and controlled model comparisons remain pending.

## Intended outcome

For the November class in Indonesia, students can follow the existing occJSDM lessons and use an AI assistant to run a supported occJSDM analysis on their own data. The preferred delivery is a Posit Cloud RStudio project with Posit Assistant and each student's OpenRouter account. A free or inexpensive model is preferred, but Doug confirmed on 5 October that a paid, more capable model is an acceptable workshop fallback. Students do not need to install R locally. This design extends the existing lessons rather than replacing them.

## What is known and what remains an assumption

Doug confirmed the class is in Indonesia, free OpenRouter accounts are acceptable, and the project should use Posit Cloud and Posit Assistant if feasible. On 4 October 2026, Doug confirmed that he successfully enabled Posit Assistant in a Posit Cloud project and shared a screenshot showing OpenRouter in the provider list. This verifies that Assistant is available in his Cloud environment and that OpenRouter is offered as a provider there. It does not yet verify that a student can connect a free individual key, that the selected free model can call MCP tools reliably, or that the Cloud environment permits the required MCP process and dependencies. Posit Assistant supports MCP server configuration. The remaining points are feasibility checks for the first pilot.

The MCP server would run in each student's Posit Cloud project container. In this sense it is local to the RStudio session, even though the project itself runs remotely in Posit Cloud. Posit Assistant supports project-level local MCP server configuration, but project settings apply only in a trusted workspace. The class setup must show students how to trust the prepared project and must check whether Posit Cloud allows the needed local subprocess. The MCP server is an R-facing interface to selected occJSDM operations, not a new statistical implementation and not arbitrary R code execution.

## Proposed architecture

Use the existing occJSDM R package as the statistical engine. Build the optional MCP component through Paper2Agent's Paper2MCP R workflow: a Python/FastMCP server invokes the original R package through Rscript. The server runs inside each student's Posit Cloud project, with isolated Python and R dependencies. Keep MCP dependencies out of the core package's required installation path. Expose a small set of named, documented tools calling existing package functions. This replaces the earlier tentative `mcptools` transport; that alternative is not part of the planned build unless the Cloud pilot identifies a concrete blocker.

Students open a prepared Posit Cloud project containing the optional MCP component, a small example dataset, and project-level MCP configuration. Each student configures their own OpenRouter key through Posit Assistant's provider setup and credential mechanism; no key is committed to the repository or project template. Posit Assistant calls the tools through MCP; the tools execute R code in that student's Cloud project and save fit artifacts there.

The assistant should receive concise validation reports and analysis summaries. The MCP tools should not transmit full datasets, arbitrary file contents, or full posterior arrays to the model by default. Model requests, tool descriptions, summaries, and any data sent to OpenRouter are subject to the provider's terms and the student's account settings; the class guide must make this clear before students use their own data.

## First class workflow

1. Check an occJSDM input object for required `info` and `OTU` structure, identify supported model structure, and report actionable issues before fitting.
2. Fit a non-spatial occJSDM model with explicit, class-safe MCMC settings and save the fit in the project.
3. Return convergence diagnostics and a compact, clearly labelled summary for discussion in Posit Assistant.

The first version targets an R object already shaped for `runOccJSDM()`, stored as an `.rds` file in the student's project. The pilot must verify this workflow with both the bundled example and a separately prepared user dataset. Additional raw-file formats, automatic covariate formula creation, maps, predictions, and broad tool catalogs are outside the first version unless the pilot shows they are essential to teach the workflow.

Do not expose spatial fitting in the initial MCP workflow. The current TODO records the spatial correction as merged, superseding this draft's original wrong-density explanation. The restriction now reflects the small classroom scope, resource budget and documented weak spatial-field recovery in some designs. Preserve the current beta limitations in interpretation. Do not expose arbitrary R evaluation or unrestricted filesystem access.

## Class and lesson integration

Doug clarified on 5 October 2026 that the eventual goal is an agent that helps each student work through the existing lesson sequence. It should use the reviewed lesson text, identify the lesson and section the student is working on, explain concepts, help with running R chunks and understanding errors, discuss the student's results, and offer questions, graduated hints and challenge problems in the student's preferred language. Let the student attempt the work and choose the pace. The four-tool quickstart pilot establishes the technical foundation for this broader teaching goal; completing it does not establish support for every lesson.

Preserve the existing lesson sequence and add a short, clearly optional AI-assisted analysis extension after the existing fitting and diagnostics lessons. Reuse existing teaching datasets and explanations where they fit. The extension should explain the MCP tool workflow, accepted input object, fit settings, how to read convergence output, and how to inspect the saved R result directly. It should show that the assistant can help operate and explain the package while the statistical estimates still come from occJSDM.

After the pilot, extend the teaching instructions and reviewed references lesson by lesson, and add only the computational tools that each lesson needs. Explanations, questions and hints belong to the conversational model and its teaching instructions; supported calculations belong to the MCP tools. Where a lesson uses computations outside the available tools, guide the student through the existing R code or use its verified saved results, clearly stating which route supplied the answer. Before claiming a lesson is supported, rehearse a student exchange covering its learning objectives, a misunderstanding or execution error, and an exercise checked against the lesson's evidence. Doug's review of the lesson text remains a prerequisite for its use as teaching material.

Provide a non-AI route through the existing lessons and R functions. It should remain possible to teach if students cannot access the selected model, encounter rate limits, or cannot use the MCP integration.

## Language and adaptive exercises

Use Bahasa Indonesia as the default language for explanations, questions, hints and feedback, with English scientific terms alongside it where helpful. Accept student responses in Indonesian, English or a mixture, and switch language on request. Maintain a reviewed bilingual glossary. Preserve R code, function names, parameter labels and data-column names exactly; translation must not alter executable arguments. Have an Indonesian scientist review representative exchanges and the terminology before the workshop.

Generate concept questions and challenge problems on demand from the lesson objectives rather than writing every question into the vignettes. Students should attempt an answer before the assistant reveals a worked answer; offer graduated hints and adapt difficulty to their responses. The first version covers conceptual questions and interpretation of saved-fit diagnostics. Numerical expected answers must come from actual R results or saved source-backed references, with the generating settings and relevant diagnostic qualifications recorded. The conversational model must not invent reference numbers or claim a short test fit establishes reliable inference. Simulation investigations and new numerical analysis operations are later additions requiring their own tested MCP tools.

These language and exercise rules belong in project-level instructions supplied by Posit Assistant to the selected model. The prepared class project will ship instruction text, the glossary and references to the saved exercise results, and setup must verify how Assistant loads them. Merely placing an instruction file beside the MCP server does not establish that the client uses it. The model generates the conversation; the MCP server runs the supported R calculations. The ordinary R teaching route remains available, and the existing lessons remain the teaching foundation.

Evaluate both tool use and Indonesian scientific explanation when selecting a model. Include questions about replication, thresholds, convergence and the distinction between baseline and site-level occupancy, as well as a challenge whose answer is checked against a saved R result. A more capable paid model is the fallback if free or inexpensive candidates do not pass the workshop trial. Record measured latency, costs and rate-limit behaviour, and arrange any funded access without distributing a shared credential in the project template.

## Feasibility gate and success criteria

Evaluate the possible cost benefit before adopting the MCP for the class: compare an inexpensive model using the repository and R directly, the same model using the MCP, and a more capable model, with actual costs measured using the same direct R route. Use matched student tasks, inputs, settings and teaching instructions; compare execution correctness, teaching quality, reliability, time and actual model charges including failed attempts and retries. Begin the direct R baselines before conversion and complete the comparison after the MCP works in Cloud. The MCP may reduce programming effort and retries, but its ability to make an inexpensive model sufficient is a hypothesis to test. If direct R already offers better workshop value, defer MCP classroom packaging and retain that route.

The first deliverable is a disposable technical pilot, not package integration. Assistant availability and the presence of OpenRouter in the provider list are verified in Doug's Posit Cloud project. The pilot must confirm a student can connect an individual OpenRouter key, select a free, inexpensive or paid model that supports tool calls, connect to an MCP server running in the project, invoke a tool, and receive the expected result. Verify a complete small non-spatial fit and diagnostic summary in the Cloud project's resource and runtime limits. Record the versions, exact setup steps, model identifier, and any limitations. If any required connection is unavailable, revise the client or hosting architecture before integrating it into the package or class project.

The class-ready workflow is successful when a student can start from the shared project, provide an eligible RDS input, get validation feedback, run the supported fit, retrieve diagnostics, and locate the saved result without editing MCP configuration by hand. The same lesson must remain teachable through the non-AI route.

## Scope boundaries

- Optional MCP support lives at root `mcp/` alongside occJSDM and does not make Python, FastMCP or an MCP library mandatory for ordinary package users. General specs and plans live in `dev/superpowers/`.
- The AI client is Posit Assistant in RStudio on Posit Cloud, contingent on the feasibility gate.
- The model provider is OpenRouter using individual API keys, with free or inexpensive models preferred and a paid model permitted for the workshop.
- The first workflow supports an occJSDM-ready RDS object and non-spatial fitting plus diagnostics.
- The MCP layer uses allowlisted tools with constrained inputs, bounded runtime settings, and project-scoped file access.
- Existing lessons remain the teaching foundation. The pilot adds one optional lesson extension and a guided setup/checklist; the eventual teaching agent will guide students through the reviewed lesson sequence in stages.
- No central hosting service, shared class API key, arbitrary R execution, spatial fitting or broad analysis assistant is part of this pilot. A paid model is an allowed workshop fallback, not a student subscription requirement.

## Risks and mitigations

Doug has verified that Posit Assistant is enabled and OpenRouter is listed in his Posit Cloud project, reducing the provider-access uncertainty. The pilot still needs to establish how student accounts connect their individual keys and whether a free model can call tools reliably. Project MCP settings may require students to trust the project, and the Cloud environment may restrict launching local subprocesses; verify both. Free OpenRouter models can change, be rate-limited, or fail to call tools consistently; pin a verified model for the lesson where possible and retain the non-AI route. Cloud CPU, memory, or runtime limits may make MCMC fits impractical; use an explicit small teaching configuration and validate runtime, and retain precomputed results for the lesson. Students may use data they should not send to a third-party model; keep raw data out of tool responses and explain what content the assistant receives. Input format mistakes can block novices; start from an RDS object in the package's native input shape and give a worked conversion example in the lesson.

## Decisions needed after the pilot

If the Posit Cloud and OpenRouter integration works, confirm the preferred and fallback models, measured resource settings, supported input structure, Indonesian teaching instructions, data handling language and class project distribution method before implementation. If it does not work, choose between RStudio Desktop plus a local MCP server, another MCP-capable client, or a non-MCP Posit Cloud workflow. No fallback is selected until the pilot identifies the actual blocker.

## References

- Posit Assistant provider setup, including OpenRouter: https://assistant.posit.co/docs/getting-started/providers/
- Posit Assistant MCP configuration, including local stdio servers and trusted project settings: https://assistant.posit.co/docs/reference/mcp-servers/
- Posit Assistant configuration reference, including project-level MCP settings: https://assistant.posit.co/docs/reference/config-file/
- Posit Cloud AI setup and its current Posit AI subscription note: https://docs.posit.co/cloud/get_started/index.html
- Paper2Agent Paper2MCP R workflow: https://github.com/jmiao24/Paper2Agent/blob/main/skills/paper2agent/paper2mcp/references/routes/r.md
- Implementation tasks: [occJSDM Paper2Agent MCP implementation plan](../plans/2026-10-05-occjsdm-paper2agent-mcp.md)
- Posit Cloud project containers and resource limits: https://docs.posit.co/cloud/guide/projects/index.html
