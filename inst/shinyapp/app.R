# taxdiv Explorer — Shiny Dashboard
# Development: feature/shiny-explorer
# ============================================================

library(shiny)
library(bslib)
library(DT)
library(shinyjs)

# ── Reusable sidebar pieces ───────────────────────────────────────────────────

sidebar_analiz <- sidebar(
  width = 320,
  bg    = "#eef0f5",

  # Load Data
  card(
    card_header("Load Data"),
    fileInput("file_upload", "Select Excel or CSV file:",
              accept      = c(".xlsx", ".xls", ".csv"),
              buttonLabel = "Browse",
              placeholder = "No file selected"),
    tags$hr(style = "margin:8px 0;"),
    tags$label("or use package dataset:"),
    selectInput("pkg_data", label = NULL,
                choices  = c("-- Select --"                      = "",
                             "anatolian_trees (correct format)"  = "anatolian_trees",
                             "gazi_comm (incorrect format example)" = "gazi_comm"),
                selected = ""),
    tags$hr(style = "margin:8px 0;"),
    downloadButton("dl_template", "Download template for data format",
                   class = "btn-outline-secondary btn-sm w-100"),
    tags$small(style = "color:#666; display:block; margin-top:4px;",
               "Arrange your data according to this template before uploading.")
  ),

  # Site Filter (visible only for multi-site data)
  uiOutput("site_filter_ui"),

  # Status box
  uiOutput("status_box"),

  # Calculation Settings
  card(
    card_header("Calculation Settings"),
    tags$label("Index groups to compute:",
               style = "font-weight:600; margin-bottom:2px;"),
    tags$small(style = "color:#666; display:block; margin-bottom:10px;",
               "Select which diversity indices to calculate.",
               "All groups are computed by default."),

    checkboxInput("idx_classical",
                  "Classical (Shannon, Simpson)", value = TRUE),
    tags$div(style = "color:#888; font-size:0.78rem; margin:-8px 0 10px 26px;",
             "Species-level diversity measures"),

    checkboxInput("idx_clarke",
                  "Clarke & Warwick (Delta, AvTD, ...)", value = TRUE),
    tags$div(style = "color:#888; font-size:0.78rem; margin:-8px 0 10px 26px;",
             "Taxonomy-weighted distinctness measures"),

    checkboxInput("idx_pto",
                  "Ozkan pTO (uTO, TO, uTO+, TO+)", value = TRUE),
    tags$div(style = "color:#888; font-size:0.78rem; margin:-8px 0 10px 26px;",
             "Deng entropy-based taxonomic originality"),

    conditionalPanel(
      condition = "input.idx_pto === true",
      tags$hr(style = "margin:6px 0 10px;"),
      numericInput("n_iter",
                   label = "Number of iterations (n_iter):",
                   value = 500,
                   min   = 1,
                   max   = 9999,
                   step  = 100),
      tags$small(style = "color:#666; display:block; margin-bottom:8px;",
                 "Default: 500. Lower values are faster, higher values give more precise estimates."),
      numericInput("seed",
                   label = "Random seed:",
                   value = NA,
                   min   = 1,
                   max   = 99999,
                   step  = 1),
      tags$small(style = "color:#666;",
                 "Leave empty for different random results each run. ",
                 "Enter a number for reproducible results (paper/review).")
    )
  ),

  # RUN button
  div(
    style = "text-align:center; margin:16px 0;",
    actionButton("btn_run", "RUN",
                 class = "btn btn-success btn-lg w-100",
                 style = "font-size:16px;")
  )
)

sidebar_grafik <- sidebar(
  width = 320,
  bg    = "#eef0f5",

  card(
    card_header("Graph Settings"),
    selectInput("graph_type",
                "Graph type:",
                choices = c(
                  "Bar Chart (Site x Index)"   = "bar",
                  "Heatmap (taxonomic distance)" = "heatmap",
                  "Taxonomic Tree (dendrogram)"  = "tree",
                  "Bubble (species contributions)" = "bubble",
                  "Radar (multi-site comparison)"  = "radar"
                ),
                selected = "bar"),

    # Bar chart settings
    conditionalPanel(
      condition = "input.graph_type == 'bar'",
      uiOutput("bar_index_ui")
    ),

    # Heatmap settings
    conditionalPanel(
      condition = "input.graph_type == 'heatmap'",
      uiOutput("heatmap_site_ui")
    ),

    # Tree settings
    conditionalPanel(
      condition = "input.graph_type == 'tree'",
      uiOutput("tree_site_ui"),
      selectInput("tree_color", "Color by:",
                  choices = c("Family", "Order", "Genus"),
                  selected = "Family")
    ),

    # Bubble settings
    conditionalPanel(
      condition = "input.graph_type == 'bubble'",
      uiOutput("bubble_site_ui"),
      selectInput("bubble_color", "Color by:",
                  choices = c("Family", "Order", "Genus"),
                  selected = "Family")
    ),

    # Radar settings (multi-site selection)
    conditionalPanel(
      condition = "input.graph_type == 'radar'",
      uiOutput("radar_sites_ui"),
      tags$small(style = "color:#666; display:block; margin-top:6px;",
                 "Tip: 2-5 sites compare best.")
    )
  ),

  card(
    card_header("Download"),
    downloadButton("dl_graph_png", "Download PNG",
                   class = "btn-outline-primary btn-sm w-100")
  )
)

# ── UI ────────────────────────────────────────────────────────────────────────
ui <- page_navbar(

  id = "main_tabs",
  theme = bs_theme(bootswatch = "flatly", base_font = font_google("Inter"),
                   font_scale = 0.85),
  header = tagList(
    tags$head(tags$title("taxdiv Explorer")),
    tags$style(HTML("
      .dataTables_wrapper { font-size: 0.82rem; }
      .card-header        { font-size: 0.88rem; padding: 6px 12px; }
      .card-body          { padding: 10px 12px; }
      label, .shiny-input-container { font-size: 0.83rem; }
      .form-control, .selectize-input { font-size: 0.83rem !important; }
      .checkbox label     { font-size: 0.82rem; }
      .nav-link           { font-weight: 500; }
      /* Selectize dropdown: roomier and not clipped by parent */
      .selectize-dropdown          { max-height: 320px !important; }
      .selectize-dropdown-content  { max-height: 300px !important;
                                     padding: 4px 0 !important; }
      .selectize-dropdown .option  { padding: 8px 12px !important;
                                     font-size: 0.85rem !important;
                                     line-height: 1.4 !important; }
      /* Let dropdowns escape ONLY where needed (sidebar cards + card
         headers). Do NOT make main-content card bodies overflow:visible,
         or tall tables (e.g. Data Preview 'All') spill out and overlap
         the panels below them. */
      .sidebar .card, .sidebar .card-body { overflow: visible !important; }
      .card-header                        { overflow: visible !important; }
    ")),
    useShinyjs()
  ),

  # ── Title (clickable: returns to Analysis tab) ────────────────────────────
  title = tags$a(
    id      = "home_link",
    href    = "#",
    class   = "action-button",
    onclick = "return false;",
    style   = "color: inherit; text-decoration: none; cursor: pointer; display: inline-flex; align-items: center;",
    title   = "Return to Analysis",
    tags$img(
      src    = "data:image/svg+xml;base64,PHN2ZyB4bWxucz0iaHR0cDovL3d3dy53My5vcmcvMjAwMC9zdmciIHZpZXdCb3g9IjAgMCAyNCAyNCIgZmlsbD0id2hpdGUiIHdpZHRoPSIzMnB4IiBoZWlnaHQ9IjMycHgiPjxwYXRoIGQ9Ik0xNyAxMmgtNVY3aC0ydjVINmwtNC41IDQuNUw2IDE2aDRsLTIgNGgxMGwtMi00aDR6Ii8+PC9zdmc+",
      height = "28px",
      style  = "margin-right:10px;"
    ),
    tags$span("taxdiv Explorer",
              style = "font-size:1.15rem; font-weight:700;")
  ),

  # ── TAB 1: Analysis ───────────────────────────────────────────────────────
  nav_panel(
    "Analysis",
    icon = icon("table"),
    page_sidebar(
      fillable = FALSE,
      sidebar  = sidebar_analiz,

      # Data Preview
      card(
        fill = FALSE,
        card_header(
          div(style = "display:flex; justify-content:space-between; align-items:center;",
              "Data Preview",
              div(style = "font-weight:normal;",
                  selectInput("preview_rows", label = NULL,
                              choices  = c("10 rows" = 10, "25 rows" = 25,
                                           "All" = 99999),
                              selected = 10,
                              width    = "120px")))
        ),
        DTOutput("preview_table")
      ),

      # Results panel
      uiOutput("results_panel")
    )
  ),

  # ── TAB 2: Graphs ─────────────────────────────────────────────────────────
  nav_panel(
    "Graphs",
    icon = icon("chart-line"),
    page_sidebar(
      fillable = FALSE,
      sidebar  = sidebar_grafik,

      uiOutput("graph_panel")
    )
  )
)

# ── SERVER ────────────────────────────────────────────────────────────────────
server <- function(input, output, session) {

  # Null-coalescing helper (used throughout)
  `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

  # Home link in header: clicking the logo/title returns to Analysis tab
  # AND resets all inputs to defaults (acts as "start fresh")
  observeEvent(input$home_link, {
    bslib::nav_select(id = "main_tabs", selected = "Analysis")

    # Clear data inputs
    shinyjs::reset("file_upload")
    updateSelectInput(session, "pkg_data", selected = "")

    # Reset index selection to defaults (all on)
    updateCheckboxInput(session, "idx_classical", value = TRUE)
    updateCheckboxInput(session, "idx_clarke",    value = TRUE)
    updateCheckboxInput(session, "idx_pto",       value = TRUE)

    # Reset pTO settings
    updateNumericInput(session, "n_iter", value = 101)
    updateNumericInput(session, "seed",   value = NA)

    # Clear results
    results(NULL)
    calc_meta(NULL)
  }, ignoreInit = TRUE)

  # "Go to Analysis" button inside graph placeholder
  observeEvent(input$go_to_analysis, {
    bslib::nav_select(id = "main_tabs", selected = "Analysis")
  }, ignoreInit = TRUE)

  # Mutual exclusion: file upload and package dataset are mutually exclusive
  # When user uploads a file, clear the package dataset dropdown
  observeEvent(input$file_upload, {
    if (!is.null(input$file_upload)) {
      updateSelectInput(session, "pkg_data", selected = "")
    }
  }, ignoreInit = TRUE)

  # When user picks a package dataset, reset the file upload
  observeEvent(input$pkg_data, {
    if (!is.null(input$pkg_data) && input$pkg_data != "") {
      shinyjs::reset("file_upload")
    }
  }, ignoreInit = TRUE)

  # Reactive: raw data
  raw_data <- reactive({
    if (!is.null(input$pkg_data) && input$pkg_data != "") {
      data(list = input$pkg_data, package = "taxdiv", envir = environment())
      return(get(input$pkg_data))
    }
    req(input$file_upload)
    ext <- tools::file_ext(input$file_upload$name)
    if (ext %in% c("xlsx", "xls")) {
      readxl::read_excel(input$file_upload$datapath)
    } else if (ext == "csv") {
      utils::read.csv(input$file_upload$datapath, stringsAsFactors = FALSE)
    }
  })

  # Detect site column
  site_col <- reactive({
    df <- raw_data()
    if (is.null(df) || !is.data.frame(df)) return(NULL)
    site_names <- c("site", "plot")
    match <- which(tolower(names(df)) %in% site_names)
    if (length(match) > 0) names(df)[match[1]] else NULL
  })

  # Site list (natural sort: Site1, Site2, ..., Site10)
  site_choices <- reactive({
    col <- site_col()
    df  <- raw_data()
    if (is.null(col)) return(NULL)
    sn       <- unique(as.character(df[[col]]))
    txt_part <- sub("(\\d+)$", "", sn)
    num_part <- suppressWarnings(as.numeric(sub("^\\D*(\\d+)$", "\\1", sn)))
    num_part[is.na(num_part)] <- 0
    sn[order(txt_part, num_part)]
  })

  # Site filter UI (only for multi-site data)
  output$site_filter_ui <- renderUI({
    choices <- site_choices()
    if (is.null(choices) || length(choices) <= 1) return(NULL)
    card(
      card_header("Site Filter"),
      div(
        style = "max-height:220px; overflow-y:auto; padding-right:4px;",
        checkboxGroupInput(
          "selected_sites",
          label    = paste0(length(choices), " sites found:"),
          choices  = choices,
          selected = choices
        )
      ),
      div(
        style = "display:flex; gap:8px; margin-top:4px;",
        actionButton("select_all",   "Select All",
                     class = "btn btn-outline-secondary btn-sm"),
        actionButton("deselect_all", "Clear",
                     class = "btn btn-outline-secondary btn-sm")
      )
    )
  })

  # Select all / clear
  observeEvent(input$select_all, {
    updateCheckboxGroupInput(session, "selected_sites", selected = site_choices())
  })
  observeEvent(input$deselect_all, {
    updateCheckboxGroupInput(session, "selected_sites", selected = character(0))
  })

  # Filtered data
  filtered_data <- reactive({
    df  <- raw_data()
    col <- site_col()
    if (is.null(df) || !is.data.frame(df) || is.null(col)) return(df)

    # If site column exists but Site Filter UI didn't render
    # (single-site dataset or filter not yet built) → use all data.
    if (is.null(input$selected_sites) || length(input$selected_sites) == 0) {
      return(df)
    }
    df[as.character(df[[col]]) %in% input$selected_sites, ]
  })

  # Preview table
  preview_data <- reactive({
    req(raw_data())
    n  <- as.integer(input$preview_rows)
    df <- raw_data()
    if (!is.data.frame(df)) {
      df <- data.frame(
        Species    = names(df),
        Abundance  = as.numeric(df),
        stringsAsFactors = FALSE
      )
    }
    head(df, if (is.na(n)) 10L else n)
  })

  output$preview_table <- renderDT(server = FALSE, {
    datatable(preview_data(),
              options = list(scrollX = TRUE, dom = "t", paging = FALSE,
                             # Cap the preview height: few rows stay short
                             # (scrollCollapse), many rows ("All") scroll
                             # inside this window instead of growing the
                             # card and overlapping the Results panel below.
                             scrollY = "420px", scrollCollapse = TRUE),
              rownames = FALSE,
              class    = "compact stripe hover")
  })

  # Status box
  output$status_box <- renderUI({
    df <- raw_data()
    if (is.null(df)) return(NULL)

    error_box <- function(title, description) {
      div(style = "background:#f8d7da; border:1px solid #f5c2c7; border-radius:8px; padding:12px 16px; margin:8px 0;",
          tags$b(paste0("Error: ", title)), tags$br(), description)
    }

    if (!is.data.frame(df)) {
      return(error_box(
        "Data is not in data frame format.",
        paste0("Loaded data is of type '", class(df), "'. ",
               "batch_analysis() requires a table with columns like Species, Genus, Family. ",
               "Try selecting the anatolian_trees dataset from the examples.")
      ))
    }

    required_cols <- c("Species", "Genus", "Family")
    missing       <- required_cols[!tolower(required_cols) %in% tolower(names(df))]

    if (length(missing) > 0) {
      return(error_box(
        paste("Missing columns:", paste(missing, collapse = ", ")),
        paste0("These columns were not found in the data: ", paste(missing, collapse = ", "), ". ",
               "Available columns: ", paste(names(df), collapse = ", "), ".")
      ))
    }

    # Abundance column check
    abd_idx <- which(tolower(names(df)) == "abundance")
    if (length(abd_idx) > 0) {
      abd_col <- names(df)[abd_idx[1]]
      if (!is.numeric(df[[abd_col]])) {
        converted <- suppressWarnings(as.numeric(df[[abd_col]]))
        bad_rows  <- which(is.na(converted) & !is.na(df[[abd_col]]))
        if (length(bad_rows) > 0) {
          return(error_box(
            paste0("'", abd_col, "' column is not numeric"),
            paste0(length(bad_rows), " row(s) contain invalid values. ",
                   "Row(s): ", paste(bad_rows, collapse = ", "), ". ",
                   "Values: ", paste(df[[abd_col]][bad_rows], collapse = ", "))
          ))
        }
      }
    }

    col      <- site_col()
    n_sites  <- if (!is.null(col)) length(unique(df[[col]])) else 1
    site_tx  <- if (n_sites > 1) paste0(", ", n_sites, " sites") else ""

    div(style = "background:#d1e7dd; border:1px solid #badbcc; border-radius:8px; padding:12px 16px; margin:8px 0;",
        tags$b("Data is valid. "),
        paste0(nrow(df), " rows, ", ncol(df), " columns", site_tx,
               " loaded. Press the button to calculate."))
  })

  # Calculation
  results    <- reactiveVal(NULL)
  calc_meta  <- reactiveVal(NULL)

  # Reset results when any input changes
  observe({
    input$file_upload
    input$pkg_data
    input$selected_sites
    input$idx_classical
    input$idx_clarke
    input$idx_pto
    input$n_iter
    input$seed
    results(NULL)
    calc_meta(NULL)
  }) |> bindEvent(
    input$file_upload, input$pkg_data, input$selected_sites,
    input$idx_classical, input$idx_clarke, input$idx_pto,
    input$n_iter, input$seed,
    ignoreInit = TRUE
  )

  observeEvent(input$btn_run, {
    req(filtered_data())
    df <- as.data.frame(filtered_data())

    if (!is.data.frame(df) || nrow(df) == 0) {
      showNotification("Please select valid data.", type = "warning")
      return()
    }

    # Build indices vector from individual checkboxes
    sel_indices <- c()
    if (isTRUE(input$idx_classical)) sel_indices <- c(sel_indices, "classical")
    if (isTRUE(input$idx_clarke))    sel_indices <- c(sel_indices, "clarke_warwick")
    if (isTRUE(input$idx_pto))       sel_indices <- c(sel_indices, "ozkan_pto")

    if (length(sel_indices) == 0) {
      showNotification("Please select at least one index group.", type = "warning")
      return()
    }

    # Auto-coerce abundance column to numeric (character -> numeric)
    abd_idx <- which(tolower(names(df)) == "abundance")
    if (length(abd_idx) > 0) {
      abd_col <- names(df)[abd_idx[1]]
      if (!is.numeric(df[[abd_col]])) {
        df[[abd_col]] <- suppressWarnings(as.numeric(df[[abd_col]]))
      }
    }

    # Determine if pTO is requested
    do_pto <- "ozkan_pto" %in% sel_indices

    # Lock button
    shinyjs::disable("btn_run")
    updateActionButton(session, "btn_run", label = "Calculating...")
    on.exit({
      shinyjs::enable("btn_run")
      updateActionButton(session, "btn_run", label = "RUN")
    })

    withProgress(message = "Calculating...", value = 0.02,
                 detail = "Preparing data", {
      used_seed <- if (do_pto) {
        if (is.na(input$seed) || is.null(input$seed)) {
          sample.int(99999L, 1L)
        } else {
          as.integer(input$seed)
        }
      } else {
        NULL
      }
      used_n_iter <- if (do_pto) as.integer(input$n_iter) else 101L

      # Per-site progress callback: bar advances after each site completes
      shiny_progress <- function(i, n, site) {
        setProgress(
          value  = 0.05 + 0.92 * (i / n),
          detail = sprintf("Site %d / %d  (%s)", i, n, site)
        )
      }

      setProgress(value = 0.05, detail = "Starting analysis...")

      t_start <- proc.time()[["elapsed"]]
      res <- tryCatch({
        taxdiv::batch_analysis(
          data        = df,
          indices     = sel_indices,
          full        = do_pto,
          n_iter      = used_n_iter,
          seed        = if (!is.null(used_seed)) used_seed else 42L,
          progress    = FALSE,
          progress_fn = shiny_progress
        )
      }, error = function(e) {
        showNotification(paste("Error:", e$message), type = "error", duration = 10)
        NULL
      })
      elapsed <- round(proc.time()[["elapsed"]] - t_start, 1)
      setProgress(value = 1, detail = "Done")
      results(res)
      calc_meta(list(
        elapsed  = elapsed,
        n_iter   = used_n_iter,
        seed     = used_seed,
        date     = format(Sys.time(), "%Y-%m-%d %H:%M"),
        indices  = sel_indices
      ))
    })
  })


  # Pick the best available index for chart/summary
  chart_index <- reactive({
    res <- results()
    if (is.null(res)) return(NULL)
    # Priority: uTO > Shannon > Delta
    candidates <- c("uTO", "Shannon", "Delta")
    for (idx in candidates) {
      if (idx %in% names(res)) return(idx)
    }
    NULL
  })

  # Results panel
  output$results_panel <- renderUI({
    req(results())
    res     <- results()
    idx_col <- chart_index()

    panels <- tagList(
      uiOutput("calc_info_bar")
    )

    # Summary boxes when a plottable index exists
    if (!is.null(idx_col)) {
      panels <- tagList(panels, uiOutput("summary_boxes"))
    }

    # Results table (always shown). Download buttons live in the card
    # header (right-aligned) so they stay visible no matter how tall the
    # table grows or which pagination page is active.
    panels <- tagList(panels,
      card(
        fill = FALSE,
        card_header(
          div(
            style = "display:flex; justify-content:space-between; align-items:center; gap:12px; flex-wrap:wrap;",
            "Results",
            div(
              style = "display:flex; gap:8px; flex-wrap:wrap; font-weight:normal;",
              downloadButton("dl_xlsx", "Download Excel (.xlsx)",
                             class = "btn-success"),
              downloadButton("dl_csv", "Download CSV (.csv)",
                             class = "btn-outline-success")
            )
          )
        ),
        DTOutput("results_table")
      )
    )

    # Bar chart moved to Graphs tab
    panels
  })

  # Calculation info bar
  output$calc_info_bar <- renderUI({
    meta <- calc_meta()
    req(meta)

    # Build index group labels
    idx_labels <- c(
      classical      = "Classical",
      clarke_warwick = "Clarke & Warwick",
      ozkan_pto      = "Ozkan pTO"
    )
    idx_text <- paste(idx_labels[meta$indices], collapse = ", ")
    has_pto  <- "ozkan_pto" %in% meta$indices

    separator <- HTML("<span style='border-left:1px solid #9fa8da; height:16px; display:inline-block; margin:0 4px;'></span>")

    info_items <- list(
      tags$span(HTML("<b>Time:</b> "), paste0(meta$elapsed, " s")),
      separator,
      tags$span(HTML("<b>Indices:</b> "), idx_text)
    )

    # Show n_iter and seed only when pTO was computed
    if (has_pto) {
      info_items <- c(info_items, list(
        separator,
        tags$span(HTML("<b>n_iter:</b> "), meta$n_iter),
        separator,
        tags$span(HTML("<b>seed:</b> "), if (!is.null(meta$seed)) meta$seed else "random")
      ))
    }

    info_items <- c(info_items, list(
      separator,
      tags$span(HTML("<b>Date:</b> "), meta$date)
    ))

    div(
      style = paste0(
        "background:#e8eaf6; border:1px solid #c5cae9; border-radius:8px; ",
        "padding:8px 16px; margin-bottom:10px; font-size:0.82rem; color:#333; ",
        "display:flex; gap:24px; flex-wrap:wrap; align-items:center;"
      ),
      info_items
    )
  })

  # Summary boxes
  output$summary_boxes <- renderUI({
    res <- results()
    req(res)

    idx_col <- chart_index()
    if (is.null(idx_col)) return(NULL)

    n_sites       <- nrow(res)
    vals          <- res[[idx_col]]
    site_col_name <- names(res)[1]

    avg  <- round(mean(vals, na.rm = TRUE), 4)
    max_ <- round(max(vals,  na.rm = TRUE), 4)
    min_ <- round(min(vals,  na.rm = TRUE), 4)

    max_site <- as.character(res[[site_col_name]][which.max(vals)])
    min_site <- as.character(res[[site_col_name]][which.min(vals)])

    box <- function(color, label, value, sub = "") {
      div(
        style = paste0(
          "background:", color, "; border-radius:10px; padding:14px 18px; ",
          "flex:1; text-align:center; min-width:140px;"
        ),
        tags$div(style = "font-size:0.78rem; color:#555; margin-bottom:4px;", label),
        tags$div(style = "font-size:1.5rem; font-weight:700;", value),
        tags$div(style = "font-size:0.75rem; color:#777; margin-top:2px;", sub)
      )
    }

    div(
      style = "display:flex; gap:12px; flex-wrap:wrap; margin-bottom:16px;",
      box("#e8f5e9", "Number of Sites", n_sites, ""),
      box("#fff9c4", paste("Mean", idx_col),    avg,     ""),
      box("#c8e6c9", paste("Highest", idx_col), max_,    max_site),
      box("#ffcdd2", paste("Lowest", idx_col),  min_,    min_site)
    )
  })

  # ── GRAPHS TAB ────────────────────────────────────────────────────────────

  # Build a numeric community vector + tax_tree from filtered_data()
  graph_inputs <- reactive({
    df <- filtered_data()
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) return(NULL)

    # Find column names (case-insensitive)
    sp_col  <- names(df)[which(tolower(names(df)) == "species")[1]]
    abd_col <- names(df)[which(tolower(names(df)) == "abundance")[1]]
    if (is.na(sp_col) || is.na(abd_col)) return(NULL)

    # Aggregate by species (sum across sites if multi-site)
    abd <- as.numeric(df[[abd_col]])
    abd[is.na(abd)] <- 0
    sp  <- as.character(df[[sp_col]])
    comm_total <- tapply(abd, sp, sum, na.rm = TRUE)

    # Build tax_tree: take first row per species
    tax_cols <- intersect(c("Genus","Family","Order","Class","Phylum","Kingdom"),
                          names(df))
    tax_tree <- df[!duplicated(df[[sp_col]]),
                   c(sp_col, tax_cols), drop = FALSE]
    names(tax_tree)[1] <- "Species"

    list(community = comm_total, tax_tree = tax_tree)
  })

  # Available index columns for bar chart
  output$bar_index_ui <- renderUI({
    res <- results()
    if (is.null(res)) return(NULL)
    site_col_name <- names(res)[1]
    idx_choices <- setdiff(names(res), c(site_col_name, "N_Species"))
    selectInput("bar_index", "Index:",
                choices  = idx_choices,
                selected = if ("uTO" %in% idx_choices) "uTO"
                           else if ("Shannon" %in% idx_choices) "Shannon"
                           else idx_choices[1])
  })

  # Helper: per-site dropdown UI generator
  site_picker_ui <- function(input_id) {
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    col <- site_col()
    if (is.null(col)) return(tags$small(style = "color:#666;",
                                        "Single-site dataset"))
    choices <- site_choices()
    selectInput(input_id, "Site:", choices = choices)
  }

  output$heatmap_site_ui <- renderUI({ site_picker_ui("heatmap_site") })
  output$tree_site_ui    <- renderUI({ site_picker_ui("tree_site") })
  output$bubble_site_ui  <- renderUI({ site_picker_ui("bubble_site") })

  # Multi-site picker for radar (checkbox list, default: first 3 sites)
  output$radar_sites_ui <- renderUI({
    df <- filtered_data()
    if (is.null(df)) return(NULL)
    col <- site_col()
    if (is.null(col)) return(tags$small(style = "color:#666;",
                                        "Single-site dataset"))
    choices <- site_choices()
    default <- head(choices, 3)
    tagList(
      tags$label(paste0("Sites to compare (", length(choices), " available):"),
                 style = "font-weight:500; margin-bottom:4px;"),
      div(
        style = paste0("max-height:220px; overflow-y:auto; padding:8px 10px; ",
                       "border:1px solid #dee2e6; border-radius:6px; ",
                       "background:white; margin-bottom:6px;"),
        checkboxGroupInput(
          "radar_sites",
          label    = NULL,
          choices  = choices,
          selected = default
        )
      ),
      div(
        style = "display:flex; gap:8px;",
        actionButton("radar_select_all", "Select All",
                     class = "btn btn-outline-secondary btn-sm"),
        actionButton("radar_clear",      "Clear",
                     class = "btn btn-outline-secondary btn-sm")
      )
    )
  })

  # Radar Select All / Clear handlers
  observeEvent(input$radar_select_all, {
    updateCheckboxGroupInput(session, "radar_sites",
                             selected = site_choices())
  })
  observeEvent(input$radar_clear, {
    updateCheckboxGroupInput(session, "radar_sites",
                             selected = character(0))
  })

  # Helper: friendly placeholder when a graph needs user input first
  placeholder_msg <- function(text) {
    ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0.5, y = 0.5,
                        label = text, size = 5, color = "#888") +
      ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) +
      ggplot2::theme_void()
  }

  # Helper: subset community + tax_tree to species of a selected site
  site_community <- function(site_name) {
    df <- filtered_data()
    col <- site_col()
    if (is.null(df) || is.null(col) || is.null(site_name)) return(NULL)
    sub <- df[as.character(df[[col]]) == site_name, ]
    sp_col  <- names(sub)[which(tolower(names(sub)) == "species")[1]]
    abd_col <- names(sub)[which(tolower(names(sub)) == "abundance")[1]]
    abd <- as.numeric(sub[[abd_col]])
    abd[is.na(abd)] <- 0
    sp  <- as.character(sub[[sp_col]])
    comm <- tapply(abd, sp, sum, na.rm = TRUE)
    comm <- comm[comm > 0]

    tax_cols <- intersect(c("Genus","Family","Order","Class","Phylum","Kingdom"),
                          names(sub))
    tt <- sub[!duplicated(sub[[sp_col]]),
              c(sp_col, tax_cols), drop = FALSE]
    names(tt)[1] <- "Species"
    tt <- tt[tt$Species %in% names(comm), , drop = FALSE]

    list(community = comm, tax_tree = tt)
  }

  # Render the graph panel based on selected type
  output$graph_panel <- renderUI({
    df <- filtered_data()
    if (is.null(df)) {
      return(div(
        style = "padding:40px; text-align:center; color:#999;",
        tags$h4("No data loaded"),
        tags$p("Please load data in the Analysis tab first.")
      ))
    }

    gt <- input$graph_type
    # All graphs require RUN first (consistency)
    if (is.null(results())) {
      return(div(
        style = "padding:40px; text-align:center; color:#666;",
        tags$h4("No analysis results yet", style = "color:#555;"),
        tags$p(
          "Switch to the Analysis tab and press RUN to compute the diversity ",
          "indices first. Then come back here to explore the graphs."
        ),
        actionButton("go_to_analysis", "Go to Analysis tab",
                     icon = icon("arrow-left"),
                     class = "btn btn-primary",
                     style = "margin-top:12px;")
      ))
    }

    # Helper for index help entries (index name in bold)
    li_idx <- function(name, text) {
      tags$li(style = "margin-bottom:6px;",
              tags$strong(name,
                          style = "font-weight:800 !important; color:#000;"),
              tags$span(HTML(paste0(" — ", text))))
    }

    bar_explanation <- tagList(
      tags$p(style = "margin:0 0 12px 0;",
        "Compares the selected index across sites. Sites are sorted from lowest to highest."),

      tags$div(style = "margin-top:6px;",
        tags$p(style = "margin:8px 0 4px 0; font-weight:600; color:#555;",
               "Classical"),
        tags$ul(style = "padding-left:22px; margin:0 0 10px 0;",
          li_idx("Shannon",
                 "information-theoretic measure of the uncertainty about which species a randomly drawn individual belongs to. Combines species richness and evenness, with value 0 for a single-species community and ln(S) for a perfectly even community of S species."),
          li_idx("Simpson",
                 "probability that two randomly drawn individuals belong to different species (the Gini-Simpson form, 1 - D). Ranges from 0 when one species dominates to 1 when the community is perfectly even, and is robust to sample size.")
        ),

        tags$p(style = "margin:8px 0 4px 0; font-weight:600; color:#555;",
               "Clarke & Warwick"),
        tags$ul(style = "padding-left:22px; margin:0 0 10px 0;",
          li_idx("Delta",
                 "abundance-weighted mean taxonomic distance across all pairs of individuals, including same-species pairs counted as distance 0. Combines abundance distribution and taxonomic spread into a single number, and decreases when one species dominates."),
          li_idx("Delta_star",
                 "abundance-weighted mean distance over different-species pairs only, with the within-species term removed. Always greater than or equal to Delta and isolates the taxonomic-distinctness signal from the effect of dominance."),
          li_idx("AvTD",
                 "average pairwise taxonomic distance computed from presence/absence only, ignoring abundance. Sample-size independent, and the most widely used Clarke & Warwick index when only species lists are available."),
          li_idx("VarTD",
                 "variance of the pairwise distances around AvTD using presence/absence only. Measures taxonomic imbalance — a high value indicates the community contains both very closely related and very distant species rather than evenly spaced taxa.")
        ),

        tags$p(style = "margin:8px 0 4px 0; font-weight:600; color:#555;",
               "Ozkan pTO (Deng-entropy based)"),
        tags$ul(style = "padding-left:22px; margin:0;",
          li_idx("uTO",
                 "unweighted Ozkan index that combines per-level Deng entropies via a multiplicative product and incorporates abundance through the slicing procedure, with all taxonomic levels weighted equally (w = 1)."),
          li_idx("TO",
                 "weighted version of uTO where each level i contributes weight w_i = i (Species = 1, Genus = 2, Family = 3, up to Kingdom = 7), so deeper taxonomic differences count more heavily than shallow ones."),
          li_idx("uTO_plus",
                 "unweighted presence-stage Ozkan index that uses only the nk = 0 slice (full species list, no abundance). The information-theoretic counterpart of AvTD."),
          li_idx("TO_plus",
                 "weighted version of uTO_plus with level weights w_i = i. Typically the largest of the four pTO values and the closest match to the Ozkan (2018) headline indicator."),
          li_idx("uTO_max / TO_max / uTO_plus_max / TO_plus_max",
                 "the same four indices computed only on informative taxonomic levels — those with Deng entropy > 0 — skipping ranks where every species sits in one group. This matches the Ozkan (2018) Excel-macro convention.")
        )
      )
    )

    descriptions <- list(
      bar     = bar_explanation,
      heatmap = "Pairwise taxonomic distance between species within the selected site. Darker cells = more distant pairs.",
      tree    = "Hierarchical taxonomic structure of the selected site as a dendrogram. Species sharing branches are taxonomically closer.",
      bubble  = "Each species in the selected site positioned by abundance and average taxonomic distance. Larger bubbles = greater diversity contribution.",
      radar   = "Compares all diversity indices across the selected sites on a single polar plot. Larger polygon area = higher overall diversity."
    )

    titles <- list(
      bar     = "Site x Index Comparison",
      heatmap = "Taxonomic Distance Heatmap",
      tree    = "Taxonomic Tree (Dendrogram)",
      bubble  = "Species Contributions (Bubble)",
      radar   = "Multi-Site Index Radar"
    )

    # All graphs use the same compact preview height
    plot_h <- 560L

    tagList(
      card(
        fill = FALSE,
        card_header(titles[[gt]]),
        plotOutput("active_graph", height = paste0(plot_h, "px"))
      ),
      card(
        fill = FALSE,
        card_header(tags$span(
          icon("circle-info"),
          tags$span("Graph Explanation", style = "margin-left:6px;")
        )),
        div(style = "padding:4px 4px; color:#444; font-size:0.88rem; line-height:1.5;",
            descriptions[[gt]])
      )
    )
  })

  # Compute a generous height for bar chart PNG download (full chart, not preview)
  bar_chart_download_height <- reactive({
    res <- results()
    if (is.null(res)) return(7)  # inches
    max(7, min(0.25 * nrow(res) + 2, 24))
  })

  # The actual graph rendering — switches based on graph_type
  active_graph_plot <- reactive({
    gt <- input$graph_type
    gi <- graph_inputs()

    switch(gt,
      bar = {
        res <- results()
        req(res, input$bar_index)
        idx_col <- input$bar_index
        site_col_name <- names(res)[1]
        vals <- res[[idx_col]]

        df_plot <- data.frame(
          site = factor(as.character(res[[site_col_name]]),
                        levels = as.character(res[[site_col_name]][order(vals)])),
          val  = vals,
          stringsAsFactors = FALSE
        )

        ggplot2::ggplot(df_plot, ggplot2::aes(x = site, y = val, fill = val)) +
          ggplot2::geom_col(width = 0.6, show.legend = FALSE) +
          ggplot2::geom_text(ggplot2::aes(label = round(val, 3)),
                             hjust = -0.15, size = 3.2) +
          ggplot2::scale_fill_gradient(low = "#ef9a9a", high = "#66bb6a") +
          ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
          ggplot2::coord_flip() +
          ggplot2::labs(x = NULL, y = idx_col) +
          ggplot2::theme_minimal(base_size = 11) +
          ggplot2::theme(
            panel.grid.major.y = ggplot2::element_blank(),
            axis.text.y        = ggplot2::element_text(size = 10),
            plot.margin        = ggplot2::margin(8, 24, 8, 8)
          )
      },

      heatmap = {
        req(gi)
        # Multi-site: wait for user to pick a site (avoid rendering huge matrix)
        if (!is.null(site_col())) {
          if (is.null(input$heatmap_site) || input$heatmap_site == "") {
            return(placeholder_msg("Select a site from the sidebar."))
          }
          sg <- site_community(input$heatmap_site)
          req(sg)
          tt <- sg$tax_tree
        } else {
          tt <- gi$tax_tree
        }
        # Adaptive label size: bigger for few species, smaller for many
        n_sp <- nrow(tt)
        lsize <- if (n_sp <= 10) 5
                 else if (n_sp <= 20) 4
                 else if (n_sp <= 35) 3
                 else 2.4
        taxdiv::plot_heatmap(tt, label_size = lsize, title = NULL)
      },

      tree = {
        req(gi)
        if (!is.null(site_col())) {
          if (is.null(input$tree_site) || input$tree_site == "") {
            return(placeholder_msg("Select a site from the sidebar."))
          }
          sg <- site_community(input$tree_site)
          req(sg)
          taxdiv::plot_taxonomic_tree(sg$tax_tree,
                                      community = sg$community,
                                      color_by   = input$tree_color %||% "Family",
                                      label_size = 3.2, title = NULL)
        } else {
          taxdiv::plot_taxonomic_tree(gi$tax_tree,
                                      community = gi$community,
                                      color_by   = input$tree_color %||% "Family",
                                      label_size = 3.2, title = NULL)
        }
      },

      bubble = {
        req(gi)
        if (!is.null(site_col())) {
          if (is.null(input$bubble_site) || input$bubble_site == "") {
            return(placeholder_msg("Select a site from the sidebar."))
          }
          sg <- site_community(input$bubble_site)
          req(sg)
          taxdiv::plot_bubble(sg$community, sg$tax_tree,
                              color_by = input$bubble_color %||% "Family",
                              title = NULL)
        } else {
          taxdiv::plot_bubble(gi$community, gi$tax_tree,
                              color_by = input$bubble_color %||% "Family",
                              title = NULL)
        }
      },

      radar = {
        res <- results()
        req(res, gi)
        df <- filtered_data()
        col <- site_col()
        req(!is.null(col))

        # Use only selected sites
        sel <- input$radar_sites
        if (is.null(sel) || length(sel) < 2) {
          return(placeholder_msg("Select at least 2 sites to compare."))
        }

        sp_col  <- names(df)[which(tolower(names(df)) == "species")[1]]
        abd_col <- names(df)[which(tolower(names(df)) == "abundance")[1]]

        comms <- lapply(sel, function(s) {
          sub <- df[as.character(df[[col]]) == s, ]
          setNames(as.numeric(sub[[abd_col]]), as.character(sub[[sp_col]]))
        })
        names(comms) <- sel

        taxdiv::plot_radar(comms, gi$tax_tree, title = NULL)
      }
    )
  })

  output$active_graph <- renderPlot({
    p <- active_graph_plot()
    print(p)
  }, height = 560)

  # PNG download (bar chart: taller to fit all sites; others: standard)
  output$dl_graph_png <- downloadHandler(
    filename = function() graph_filename(),
    content = function(file) {
      p <- active_graph_plot()
      h_in <- if (input$graph_type == "bar") bar_chart_download_height() else 7
      ggplot2::ggsave(file, plot = p, width = 10, height = h_in,
                      dpi = 150, bg = "white")
    }
  )

  output$results_table <- renderDT({
    req(results())
    datatable(results(),
              options = list(scrollX = TRUE, pageLength = 10,
                             lengthMenu = c(10, 25, 50, 100)),
              rownames = FALSE,
              class    = "compact stripe hover")
  })

  # Helper: shared base name of the current input (file upload or pkg dataset)
  input_basename <- function() {
    if (!is.null(input$file_upload) && !is.null(input$file_upload$name)) {
      tools::file_path_sans_ext(input$file_upload$name)
    } else if (!is.null(input$pkg_data) && input$pkg_data != "") {
      input$pkg_data
    } else {
      paste0("taxdiv_", Sys.Date())
    }
  }

  # "RESULT_<input-name>.<ext>"
  result_filename <- function(ext) {
    paste0("RESULT_", input_basename(), ".", ext)
  }

  # "GRAPH_<input-name>_<type>_<context>.png"
  graph_filename <- function() {
    gt   <- input$graph_type
    base <- input_basename()
    context <- switch(gt,
      bar     = input$bar_index %||% "index",
      heatmap = input$heatmap_site %||% "single",
      tree    = input$tree_site    %||% "single",
      bubble  = input$bubble_site  %||% "single",
      radar   = {
        n <- length(input$radar_sites %||% character(0))
        if (n == 0) "nosites" else paste0(n, "sites")
      },
      "graph"
    )
    paste0("GRAPH_", base, "_", gt, "_", context, ".png")
  }

  output$dl_xlsx <- downloadHandler(
    filename = function() result_filename("xlsx"),
    content  = function(file) openxlsx::write.xlsx(results(), file)
  )

  output$dl_csv <- downloadHandler(
    filename = function() result_filename("csv"),
    content  = function(file) utils::write.csv(results(), file, row.names = FALSE)
  )

  output$dl_template <- downloadHandler(
    filename = function() "taxdiv_data_template.xlsx",
    content  = function(file) {
      # -- ENTER_DATA: numeric coding example (3 sites, Westhoff-Maarel 1-9) --
      enter_data <- data.frame(
        Site      = c(rep("Site_1", 6), rep("Site_2", 5), rep("Site_3", 4)),
        Species   = c(1,2,3,4,5,6, 1,7,3,8,9, 2,4,10,5),
        Genus     = c(1,2,3,4,5,6, 1,7,3,8,9, 2,4,10,5),
        Family    = c(1,2,2,1,3,4, 1,1,2,5,6, 2,1,7,3),
        Order     = c(1,2,2,1,2,1, 1,1,2,3,4, 2,1,5,2),
        Class     = c(1,2,2,1,2,1, 1,1,2,2,2, 2,1,2,2),
        Phylum    = c(1,2,2,1,2,1, 1,1,2,2,2, 2,1,2,2),
        Kingdom   = rep(1, 15),
        Abundance = c(7,3,2,5,3,1, 4,5,5,3,4, 8,2,1,1),
        stringsAsFactors = FALSE
      )

      # -- EXAMPLE_text: text taxonomy names (Westhoff-Maarel 1-9) --
      example_text <- data.frame(
        Site = c(rep("Site_1", 6), rep("Site_2", 5), rep("Site_3", 4)),
        Species = c("Pinus nigra","Quercus cerris","Fagus orientalis",
                     "Cedrus libani","Carpinus betulus","Juniperus oxycedrus",
                     "Pinus nigra","Abies nordmanniana","Fagus orientalis",
                     "Tilia tomentosa","Acer platanoides",
                     "Quercus cerris","Cedrus libani","Sorbus torminalis","Carpinus betulus"),
        Genus = c("Pinus","Quercus","Fagus","Cedrus","Carpinus","Juniperus",
                   "Pinus","Abies","Fagus","Tilia","Acer",
                   "Quercus","Cedrus","Sorbus","Carpinus"),
        Family = c("Pinaceae","Fagaceae","Fagaceae","Pinaceae","Betulaceae","Cupressaceae",
                    "Pinaceae","Pinaceae","Fagaceae","Malvaceae","Sapindaceae",
                    "Fagaceae","Pinaceae","Rosaceae","Betulaceae"),
        Order = c("Pinales","Fagales","Fagales","Pinales","Fagales","Pinales",
                   "Pinales","Pinales","Fagales","Malvales","Sapindales",
                   "Fagales","Pinales","Rosales","Fagales"),
        Class = c("Pinopsida","Magnoliopsida","Magnoliopsida","Pinopsida","Magnoliopsida","Pinopsida",
                   "Pinopsida","Pinopsida","Magnoliopsida","Magnoliopsida","Magnoliopsida",
                   "Magnoliopsida","Pinopsida","Magnoliopsida","Magnoliopsida"),
        Phylum = c("Pinophyta","Magnoliophyta","Magnoliophyta","Pinophyta","Magnoliophyta","Pinophyta",
                    "Pinophyta","Pinophyta","Magnoliophyta","Magnoliophyta","Magnoliophyta",
                    "Magnoliophyta","Pinophyta","Magnoliophyta","Magnoliophyta"),
        Kingdom = rep("Plantae", 15),
        Abundance = c(7,3,2,5,3,1, 4,5,5,3,4, 8,2,1,1),
        stringsAsFactors = FALSE
      )

      # -- EXAMPLE_numeric: numeric taxonomy codes (Westhoff-Maarel 1-9) --
      example_numeric <- data.frame(
        Site      = c(rep("Site_1", 6), rep("Site_2", 5), rep("Site_3", 4)),
        Species   = c(1,2,3,4,5,6, 1,7,3,8,9, 2,4,10,5),
        Genus     = c(1,2,3,4,5,6, 1,7,3,8,9, 2,4,10,5),
        Family    = c(1,2,2,1,3,4, 1,1,2,5,6, 2,1,7,3),
        Order     = c(1,2,2,1,2,1, 1,1,2,3,4, 2,1,5,2),
        Class     = c(1,2,2,1,2,1, 1,1,2,2,2, 2,1,2,2),
        Phylum    = c(1,2,2,1,2,1, 1,1,2,2,2, 2,1,2,2),
        Kingdom   = rep(1, 15),
        Abundance = c(7,3,2,5,3,1, 4,5,5,3,4, 8,2,1,1),
        stringsAsFactors = FALSE
      )

      # -- INSTRUCTIONS --
      na <- "NA"
      instructions <- data.frame(
        Topic = c("HOW TO USE","Step 1:","Step 2:","Step 3:",na,
                   "COLUMN DESCRIPTIONS","Site (optional):","Species (required):",
                   "Genus ~ Kingdom:","Abundance (required):",na,
                   "NUMERIC CODING RULE","Rule:","Example:",na,
                   "R CODE",na,na,na,"SINGLE SITE",na),
        Description = c(na,
          "Enter your data under the headers in the ENTER_DATA sheet",
          "Save the file as .xlsx",
          "Run the R code below",na,na,
          "Site/plot name (leave blank for single-site analysis)",
          "Species name or numeric code (each row should be one species)",
          "Taxonomic rank — text names or numeric codes can be used",
          "Number of individuals (must be a numeric value)",na,
          "You can use numbers instead of names for all taxonomic columns (including Species)",
          "Species in the SAME group must have the SAME number",
          "If Pinus and Cedrus share the same Family (Pinaceae), both should have Family=1",na,
          "library(taxdiv)","library(readxl)",
          'data <- as.data.frame(read_excel("your_file.xlsx"))',
          "result <- batch_analysis(data)",
          "Leave the Site column blank or remove it entirely",
          "batch_analysis() will automatically treat it as a single site"),
        stringsAsFactors = FALSE
      )

      # -- Styles --
      hdr_data  <- openxlsx::createStyle(fontSize = 12, textDecoration = "bold",
                                          halign = "center", valign = "center")
      cell_data <- openxlsx::createStyle(fontSize = 11,
                                          halign = "center", valign = "center")
      hdr_inst  <- openxlsx::createStyle(fontSize = 12, textDecoration = "bold",
                                          halign = "left", valign = "center")
      sect_inst <- openxlsx::createStyle(fontSize = 11, textDecoration = "bold",
                                          halign = "left", valign = "center")
      cell_inst <- openxlsx::createStyle(fontSize = 11,
                                          halign = "left", valign = "center")

      wb <- openxlsx::createWorkbook()

      # -- Write data sheets with styles --
      for (sn in c("ENTER_DATA", "EXAMPLE_text", "EXAMPLE_numeric")) {
        df <- switch(sn,
          ENTER_DATA      = enter_data,
          EXAMPLE_text    = example_text,
          EXAMPLE_numeric = example_numeric
        )
        openxlsx::addWorksheet(wb, sn)
        openxlsx::writeData(wb, sn, df)
        openxlsx::addStyle(wb, sn, hdr_data, rows = 1, cols = 1:9, gridExpand = TRUE)
        openxlsx::addStyle(wb, sn, cell_data, rows = 2:(nrow(df) + 1),
                           cols = 1:9, gridExpand = TRUE)
        openxlsx::setColWidths(wb, sn, cols = 1:9, widths = 15.12)
      }

      # -- Write INSTRUCTIONS with styles --
      openxlsx::addWorksheet(wb, "INSTRUCTIONS")
      openxlsx::writeData(wb, "INSTRUCTIONS", instructions)
      openxlsx::addStyle(wb, "INSTRUCTIONS", hdr_inst, rows = 1, cols = 1:2,
                         gridExpand = TRUE)
      openxlsx::addStyle(wb, "INSTRUCTIONS", cell_inst,
                         rows = 2:(nrow(instructions) + 1), cols = 1:2,
                         gridExpand = TRUE)
      # Bold section headers: HOW TO USE, COLUMN DESCRIPTIONS,
      #   NUMERIC CODING RULE, R CODE, SINGLE SITE
      bold_rows <- c(2, 7, 13, 17, 21)
      openxlsx::addStyle(wb, "INSTRUCTIONS", sect_inst,
                         rows = bold_rows, cols = 1, gridExpand = TRUE)
      openxlsx::setColWidths(wb, "INSTRUCTIONS", cols = 1:2,
                             widths = c(18.62, 64.95))

      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
    }
  )

}

# ── RUN ───────────────────────────────────────────────────────────────────────
shinyApp(ui, server)
