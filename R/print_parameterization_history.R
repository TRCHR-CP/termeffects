# ============================================================
# Parameterization history printer
# Written by Lukas Bolsinger for p-ggeffects; copied unchanged except that
# box-drawing characters are written as \u escapes (required by CRAN).
# ============================================================

.print_parameterization_history <- function(history) {

  dag <- .build_parameterization_dag(history)

  .print_parameterization_dag(dag)

  invisible(NULL)
}


# ============================================================
# Build DAG
# ============================================================

.build_parameterization_dag <- function(history) {

  nodes <- list()
  edges <- list()

  for (stage in seq_along(history)) {

    parameterizations <- history[[stage]]

    for (nm in names(parameterizations)) {

      fun <- parameterizations[[nm]]

      # Variable introduced at this stage
      nodes[[nm]] <- list(
        name = nm,
        stage = stage
      )

      # Variables used by this parameterization
      parameters <- names(formals(fun))

      expression <- paste(
        deparse(.unwrap_body(body(fun))),
        collapse = ""
      )

      for (parameter in parameters) {

        edges[[length(edges) + 1L]] <- list(
          from = nm,
          to = parameter,
          stage = stage,
          expression = expression
        )
      }
    }
  }

  # Variables which are never parameterized are inputs.
  parameterized <- names(nodes)

  inputs <- unique(
    vapply(
      edges,
      function(edge) edge$to,
      character(1)
    )
  )

  inputs <- setdiff(
    inputs,
    parameterized
  )

  list(
    nodes = nodes,
    edges = edges,
    inputs = inputs,
    n_stages = length(history)
  )
}


# ============================================================
# Graph utilities
# ============================================================

.dag_edges_from <- function(dag, name) {

  Filter(
    function(edge) identical(edge$from, name),
    dag$edges
  )
}


.dag_edges_to <- function(dag, name) {

  Filter(
    function(edge) identical(edge$to, name),
    dag$edges
  )
}


.dag_parents <- function(dag, name) {

  edges <- .dag_edges_to(
    dag,
    name
  )

  unique(
    vapply(
      edges,
      function(edge) edge$from,
      character(1)
    )
  )
}


.dag_children <- function(dag, name) {

  edges <- .dag_edges_from(
    dag,
    name
  )

  unique(
    vapply(
      edges,
      function(edge) edge$to,
      character(1)
    )
  )
}


# ============================================================
# Determine row ordering
# ============================================================

.dag_rows <- function(dag) {

  nodes <- names(dag$nodes)

  # Start with Stage 1 nodes in their history order.
  stage1 <- nodes[
    vapply(
      dag$nodes,
      function(x) identical(x$stage, 1L),
      logical(1)
    )
  ]

  rows <- setNames(
    seq_along(stage1),
    stage1
  )

  # Add later-stage nodes.
  #
  # A node depending on several previous nodes is placed
  # roughly between its parents.
  for (stage in seq.int(
    2L,
    dag$n_stages
  )) {

    stage_nodes <- nodes[
      vapply(
        dag$nodes,
        function(x) identical(x$stage, stage),
        logical(1)
      )
    ]

    for (nm in stage_nodes) {

      parents <- .dag_parents(
        dag,
        nm
      )

      parent_rows <- rows[parents]

      parent_rows <- parent_rows[
        !is.na(parent_rows)
      ]

      if (length(parent_rows)) {
        desired <- round(
          mean(parent_rows)
        )
      } else {
        desired <- length(rows) + 1L
      }

      # Find nearest unused row.
      used <- unname(rows)

      while (desired %in% used) {
        desired <- desired + 1L
      }

      rows[[nm]] <- desired
    }
  }

  # Ensure rows are consecutive.
  order_rows <- order(
    rows,
    names(rows)
  )

  new_rows <- seq_along(order_rows)

  names(new_rows) <- names(rows)[order_rows]

  new_rows
}


# ============================================================
# Formatting
# ============================================================

.edge_label <- function(edge) {

  paste0(
    edge$from,
    " [",
    edge$expression,
    "]"
  )
}


.stage_width <- function(dag, stage) {

  nodes <- names(dag$nodes)[
    vapply(
      dag$nodes,
      function(node) {
        identical(
          node$stage,
          stage
        )
      },
      logical(1)
    )
  ]

  if (!length(nodes)) {
    return(1L)
  }

  # Width required for the longest expression in this stage.
  edges <- Filter(
    function(edge) {
      identical(
        edge$stage,
        stage
      )
    },
    dag$edges
  )

  node_width <- max(
    nchar(nodes)
  )

  edge_width <- if (length(edges)) {
    max(
      nchar(
        vapply(
          edges,
          .edge_label,
          character(1)
        )
      )
    )
  } else {
    0L
  }

  max(
    node_width,
    edge_width
  ) + 4L
}


# ============================================================
# Render one row
# ============================================================

.render_stage_row <- function(
    dag,
    rows,
    row,
    stage,
    widths) {

  nodes <- names(dag$nodes)[
    vapply(
      dag$nodes,
      function(node) {
        identical(
          node$stage,
          stage
        ) &&
          identical(
            rows[[node$name]],
            row
          )
      },
      logical(1)
    )
  ]

  edges <- Filter(
    function(edge) {
      identical(
        edge$stage,
        stage
      ) &&
        identical(
          rows[[edge$from]],
          row
        )
    },
    dag$edges
  )

  if (!length(nodes) && !length(edges)) {
    return("")
  }

  if (length(nodes)) {
    node <- nodes[[1L]]
  } else {
    node <- NULL
  }

  # A parameterization is displayed as:
  #
  # x [cos(t)]
  #
  # rather than separately displaying the node and edge.
  if (!is.null(node)) {

    text <- node

    edge <- Filter(
      function(edge) {
        identical(
          edge$from,
          node
        )
      },
      edges
    )

    if (length(edge)) {

      text <- paste0(
        node,
        " [",
        edge[[1L]]$expression,
        "]"
      )
    }

  } else {

    text <- ""
  }

  text
}


# ============================================================
# Main renderer
# ============================================================

.print_parameterization_dag <- function(dag) {

  rows <- .dag_rows(dag)

  widths <- vapply(
    seq_len(dag$n_stages),
    function(stage) {
      .stage_width(
        dag,
        stage
      )
    },
    integer(1)
  )

  nrows <- max(
    1L,
    length(rows)
  )

  # ----------------------------------------------------------
  # Header
  # ----------------------------------------------------------

  header <- paste(
    vapply(
      seq_len(dag$n_stages),
      function(stage) {
        sprintf(
          "%-*s",
          widths[[stage]],
          paste0(
            "Stage ",
            stage
          )
        )
      },
      character(1)
    ),
    collapse = "   "
  )

  separator <- paste(
    vapply(
      widths,
      function(width) {
        strrep(
          "-",
          width
        )
      },
      character(1)
    ),
    collapse = "   "
  )

  cat(
    "Parameterization history:\n\n",
    header,
    "\n",
    separator,
    "\n\n",
    sep = ""
  )

  # ----------------------------------------------------------
  # Body
  # ----------------------------------------------------------

  for (row in seq_len(nrows)) {

    pieces <- character(
      dag$n_stages
    )

    for (stage in seq_len(dag$n_stages)) {

      pieces[[stage]] <- .render_stage_row(
        dag,
        rows,
        row,
        stage,
        widths
      )
    }

    # --------------------------------------------------------
    # Add arrows between stages.
    # --------------------------------------------------------

    for (stage in seq_len(dag$n_stages - 1L)) {

      # Find whether something on this row has a dependency
      # that lands in the next stage.
      current <- names(dag$nodes)[
        vapply(
          dag$nodes,
          function(node) {
            identical(
              node$stage,
              stage
            ) &&
              identical(
                rows[[node$name]],
                row
              )
          },
          logical(1)
        )
      ]

      if (length(current)) {

        children <- unique(
          unlist(
            lapply(
              current,
              function(x) {
                .dag_children(
                  dag,
                  x
                )
              }
            )
          )
        )

        children <- children[
          children %in% names(dag$nodes)
        ]

        if (length(children)) {

          target_stage <- vapply(
            dag$nodes[children],
            function(x) {
              x$stage
            },
            integer(1)
          )

          if (any(
            target_stage == stage + 1L
          )) {

            pieces[[stage]] <- paste0(
              pieces[[stage]],
              " <\u2500\u2500\u2500\u2500\u2500\u2510"
            )
          }
        }
      }
    }

    # --------------------------------------------------------
    # Terminal inputs
    # --------------------------------------------------------

    inputs_here <- character()

    for (edge in dag$edges) {

      if (!identical(
        rows[[edge$from]],
        row
      )) {
        next
      }

      if (edge$to %in% dag$inputs) {
        inputs_here <- c(
          inputs_here,
          edge$to
        )
      }
    }

    inputs_here <- unique(
      inputs_here
    )
    if (length(inputs_here)) {

      pieces[[dag$n_stages]] <- paste0(
        pieces[[dag$n_stages]],
        " <\u2500\u2500 ",
        paste(
          inputs_here,
          collapse = ", "
        )
      )
    }
    # if (length(inputs_here)) {
    #
    #   pieces <- c(
    #     pieces,
    #     paste0(
    #       "<---- ",
    #       paste(
    #         inputs_here,
    #         collapse = ", "
    #       )
    #     )
    #   )
    # }

    # --------------------------------------------------------
    # Pad stages and print.
    # --------------------------------------------------------

    line <- paste(
      vapply(
        seq_along(pieces),
        function(i) {

          if (i <= length(widths)) {

            sprintf(
              "%-*s",
              widths[[i]],
              pieces[[i]]
            )

          } else {
            pieces[[i]]
          }
        },
        character(1)
      ),
      collapse = "   "
    )

    cat(
      sub(
        "\\s+$",
        "",
        line
      ),
      "\n",
      sep = ""
    )
  }

  invisible(NULL)
}
