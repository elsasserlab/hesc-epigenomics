#' Code adapted from Simon's embed last plot snippets

#' Convert table to text, as it would be written in a file
#'
#' @param t data.frame or anything that can be written with write.table
#'
#' @return A string
table_to_string <- function(t) {
  tmp <- tempfile("last_plot_table", fileext = "csv")
  write.table(t, tmp, quote = F, col.names = T, row.names = F, sep ="\t")
  table_string <- paste0(readLines(tmp), collapse="\n")
  unlink(tmp)
  table_string
}

#' Create encoded Base64 datastream to embed in html
#'
#' @param x data.frame or table-like object.
#'
#' @return A Base64 encoded string
encode_data <- function(x) {
  encoded_table <- openssl::base64_encode(table_to_string(x))
  sprintf('data:text/csv;base64,%s', encoded_table)
}

#' Embed last plot data in html

#' @return A html paragraph with embedded data
embed_last_plot_data <- function(name = "data.tsv", include_preview = F) {
  last_plot_table <- ggplot_build(last_plot())$plot$data
  link <- paste0("<a download='", name ,
                 "' href=", encode_data(last_plot_table), ">",
                 "download plot data","</a>")

  if (include_preview == TRUE) {
    last_plot_preview <- paste(kable(head(last_plot_table),
                                     format = "html"), collapse = "\n")
    paste("<p>Plot data: ",
          link,
          " - data table preview: </p>",
          last_plot_preview,
          collapse = "\n")
  } else {
    link
  }
}

embed_df <- function(df, name = "data.tsv") {
  paste0("<a download='",
         name,
         "' href=",
         encode_data(df),
         ">",
         "download plot data",
         "</a>")
}
