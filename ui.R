library(shiny)
#library(shinyjs)

## Note: This objects are not reloaded when the 'Reload App' button is clicked inside RStudio.
subtitlesStyle <- "font-weight: bold; color: #000000;"

ltys <- list(Solid = 1, Dashed = "33", Dotted = "12", `Dot-dash` = "2393", `Long dash` = 5, `Two dash` = 6)

# Prepairing help
tabs <- c("Polygons", "Labels", "Guides", "Legend", "Groups", "Title", "Net", "Wheel", "Style", "Tools", "Exporting")
helpTexts <- lapply(tabs, function(tab) readLines(paste0("www/help", tab, ".txt"), warn = FALSE))

## 

seqInputRow <- function(inputId, label, value = "", width) {
  div(style = "margin-bottom: 0;",
      tags$label(label, `for` = inputId), 
      tags$input(id = inputId, type = "text", value = value, width = width,
                 class="form-control shinyjs-resettable shiny-bound-input"))
}


shinyUI(fixedPage(
  tags$head(
    HTML("<script async src=\"https://www.googletagmanager.com/gtag/js?id=UA-60942909-3\"></script>"),
    includeScript("google_analytics.js"),
    tags$style(HTML("
                    hr {
                          margin: 0;
                    }
                    html {
                          overflow-y: scroll; 
                    }
                    .nav li a {
                          padding: 10px 10px;
                    }
                    ")
    ) #overflow-y forces vertical scrollbar, avoids movement when switching from Net to Wheel
  ),
  shinyjs::useShinyjs(),
  fixedPage(
    
  ),
  titlePanel("NetWheels: Peptides Helical Wheel and Net projections maker"),
  fixedPage(
    fixedRow(
      id = "side-panel",
      style = "height: 760px;",
      tags$br(),
      column(4,
             helpText("Type or paste the peptide sequence on the \"Sequence\" field. Navigate the through the tabs below to adjust several configurations of the projection. Check the \"Help\" tag above the image for detailed instructions."),
             #
             seqInputRow(inputId = "seq", label = "Sequence", value = "DLISGLGQRNVXKVLTETGLP", #QGAMNKALELFRKDI
                         width = "100%"),#, style = "margin-bottom: 0px;"),
             #textInput(inputId = "seq", label = "Sequence", value = "DLISGLGQRNVXKVLTETGLP", width = "100%"),
             uiOutput("resCount"),
             uiOutput("seqMono"),
             #uiOutput("resNumber"),
             #            numericInput(inputId = "wheelsize", label = "Wheel Size (Disabled)", 
             #                         value = 10, min = 1, max = 20, step = 0.5, width = "100%"),
             #            fixedRow(
             #              column(12, # Deixando a column caso apareça algum botao para colocar do lado
             #                     actionButton("reset_input", "Reset ALL inputs to default values")
             #              )
             #           ),
             tags$div(id = "settings",
                      tabsetPanel(
                        tabPanel("Polygons",
                                 tabsetPanel(
                                   tabPanel("Shape",
                                            fixedRow(
                                              column(6,
                                                     selectInput(inputId = "shp1", label = "Polar / Basic",
                                                                 choices = c("Circle", "Square", "Triangle", "Diamond", "Hexagon"), selected = "Circle")
                                              ),
                                              column(6,
                                                     selectInput(inputId = "shp2", label = "Polar / Acidic",
                                                                 choices = c("Circle", "Square", "Triangle", "Diamond", "Hexagon"), selected = "Circle")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     selectInput(inputId = "shp3", label = "Polar / Uncharged",
                                                                 choices = c("Circle", "Square", "Triangle", "Diamond", "Hexagon"), selected = "Circle")
                                              ),
                                              column(6,
                                                     selectInput(inputId = "shp4", label = "Nonpolar",
                                                                 choices = c("Circle", "Square", "Triangle", "Diamond", "Hexagon"), selected = "Circle")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     selectInput(inputId = "shp5", label = "Unknown Residue",
                                                                 choices = c("Circle", "Square", "Triangle", "Diamond", "Hexagon"), selected = "Circle")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     numericInput(inputId = "circprop", label = "First/Last Ratio", 
                                                                  value = 1, min = 0.1, max = 5, step = 0.05)
                                              ),
                                              column(6,
                                                     numericInput(inputId = "circsize", label = "Size", 
                                                                  value = 0.3, min = 0.1, max = 1, step = 0.01)
                                              )
                                            )
                                   ),
                                   tabPanel("Colors",
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "col1", label = "Polar / Basic", value = "red")
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "col2", label = "Polar / Acidic", value = "blue")
                                              )),
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "col3", label = "Polar / Uncharged", value = "green")
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "col4", label = "Nonpolar", value = "yellow")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "col5", label = "Unknown Residue",
                                                                          value = "white")
                                              )
                                            )
                                   ),
                                   tabPanel("Pattern",
                                            helpText("These are available only for circles. There's no reason to mix grids and different polygons."),
                                            helpText("n = none, v = vertical, h = horizontal, d = diagonals"),
                                            fixedRow(
                                              column(3,
                                                     selectInput(inputId = "fill1", label = "Basic",
                                                                 choices = c("n", "v", "h", "d/", "d\\"), selected = "n")
                                              ),
                                              column(3,
                                                     numericInput(inputId = "nFills1", label = "Lines",
                                                                  value = 3, min = 1, max = 21, step = 2)
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "fillCol1", label = "Line Color", value = "gray30")
                                              )
                                            ),
                                            fixedRow(
                                              column(3,
                                                     selectInput(inputId = "fill2", label = "Acidic",
                                                                 choices = c("n", "v", "h", "d/", "d\\"), selected = "n")
                                              ),
                                              column(3,
                                                     numericInput(inputId = "nFills2", label = "Lines",
                                                                  value = 3, min = 1, max = 21, step = 2)
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "fillCol2", label = "Line Color", value = "gray30")
                                              )
                                            ),
                                            fixedRow(
                                              column(3,
                                                     selectInput(inputId = "fill3", label = "Uncharged",
                                                                 choices = c("n", "v", "h", "d/", "d\\"), selected = "n")
                                              ),
                                              column(3,
                                                     numericInput(inputId = "nFills3", label = "Lines",
                                                                  value = 3, min = 1, max = 21, step = 2)
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "fillCol3", label = "Line Color", value = "gray30")
                                              )
                                            ),
                                            fixedRow(
                                              column(3,
                                                     selectInput(inputId = "fill4", label = "Nonpolar",
                                                                 choices = c("n", "v", "h", "d/", "d\\"), selected = "n")
                                              ),
                                              column(3,
                                                     numericInput(inputId = "nFills4", label = "Lines",
                                                                  value = 3, min = 1, max = 21, step = 2)
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "fillCol4", label = "Line Color", value = "gray30")
                                              )
                                            ),
                                            fixedRow(
                                              column(3,
                                                     selectInput(inputId = "fill5", label = "Unknown",
                                                                 choices = c("n", "v", "h", "d/", "d\\"), selected = "n")
                                              ),
                                              column(3,
                                                     numericInput(inputId = "nFills5", label = "Lines",
                                                                  value = 4, min = 1, max = 21, step = 2)
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "fillCol5", label = "Line Color", value = "gray30")
                                              )
                                            )
                                   ),
                                   tabPanel("Borders",
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "circBorder1", label = "Polar / Basic",
                                                                 choices = c("Yes", "No"), selected = "No")
                                              ),
                                              column(4,
                                                     colourpicker:: colourInput(inputId = "circBorderCol1",
                                                                          label = "Color", value = "#000000")
                                              ),
                                              column(4,
                                                     numericInput(inputId = "circBorderWd1", label = "Width", 
                                                                  value = 1, min = 0, max = 10, step = 1)
                                              )
                                            ),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "circBorder2", label = "Polar / Acidic",
                                                                 choices = c("Yes", "No"), selected = "No")
                                              ),
                                              column(4,
                                                     colourpicker:: colourInput(inputId = "circBorderCol2",
                                                                          label = "Color", value = "#000000")
                                              ),
                                              column(4,
                                                     numericInput(inputId = "circBorderWd2", label = "Width", 
                                                                  value = 1, min = 0, max = 10, step = 1)
                                              )
                                            ),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "circBorder3", label = "Polar / Uncharged",
                                                                 choices = c("Yes", "No"), selected = "No")
                                              ),
                                              column(4,
                                                     colourpicker:: colourInput(inputId = "circBorderCol3",
                                                                          label = "Color", value = "#000000")
                                              ),
                                              column(4,
                                                     numericInput(inputId = "circBorderWd3", label = "Width", 
                                                                  value = 1, min = 0, max = 10, step = 1)
                                              )
                                            ),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "circBorder4", label = "Nonpolar ",
                                                                 choices = c("Yes", "No"), selected = "No")
                                              ),
                                              column(4,
                                                     colourpicker:: colourInput(inputId = "circBorderCol4",
                                                                          label = "Color", value = "#000000")
                                              ),
                                              column(4,
                                                     numericInput(inputId = "circBorderWd4", label = "Width", 
                                                                  value = 1, min = 0, max = 10, step = 1)
                                              )
                                            ),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "circBorder5", label = "Group 5",
                                                                 choices = c("Yes", "No"), selected = "Yes")
                                              ),
                                              column(4,
                                                     colourpicker:: colourInput(inputId = "circBorderCol5",
                                                                          label = "Color", value = "#000000")
                                              ),
                                              column(4,
                                                     numericInput(inputId = "circBorderWd5", label = "Width", 
                                                                  value = 1, min = 0, max = 10, step = 1)
                                              )
                                            )
                                   )
                                 )
                        ),
                        tabPanel("Labels",
                                 tabsetPanel(
                                   tabPanel("Amino Acids",
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "labType", label = "Code Style", 
                                                                 choices = c(`1-Letter` = 1, `3-Letter` = 2,
                                                                             LetPos = 3, Position = 4, None = 0),
                                                                 selected = 1)
                                              ),
                                              column(3, 
                                                     numericInput(inputId = "labCex", label = "Font Size", 
                                                                  value = 12, min = 1, max = 30, step = 1)
                                              ),
                                              column(5, #1=plain, 2=bold, 3=italic, 4=bold italic
                                                     selectInput(inputId = "labFont", label = "Font Style", 
                                                                 choices = c(Normal = 1, Bold = 2, Italic = 3, `Bold Italic` = 4), 
                                                                 selected = 2)
                                              )
                                            ),
                                            fixedRow(
                                              column(6, 
                                                     numericInput(inputId = "labOffY", label = "Vertical Offset", 
                                                                  value = 0, min = -3, max = 2, step = 0.1)
                                              ),
                                              column(6, 
                                                     numericInput(inputId = "labOffX", label = "Horizontal Offset", 
                                                                  value = 0, min = -2, max = 2, step = 0.1)
                                              )#,
                                              #                                               column(4, 
                                              #                                                      numericInput(inputId = "labOffPos", label = "Position Offset", 
                                              #                                                                   value = 0, min = 0, max = 10000, step = 1)
                                              #                                              )
                                            ),
                                            tags$hr(),
                                            helpText("Font Color", style = subtitlesStyle),
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "labCol1", label = "Polar / Basic", value = "black")
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "labCol2", label = "Polar / Acidic", value = "white")
                                              )),
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "labCol3", label = "Polar / Uncharged", value = "black")
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "labCol4", label = "Nonpolar", value = "black")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "labCol5", label = "Unknown Residue", value = "black")
                                              )
                                            )
                                   ),
                                   tabPanel("Position",
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "numShow", label = "Show on...", 
                                                                 choices = c("Wheel", "Net", "Both", "None"), selected = "Wheel")
                                              ),
                                              column(3, 
                                                     numericInput(inputId = "numCex", label = "Font Size", 
                                                                  value = 12, min = 1, max = 30, step = 1)
                                              ),
                                              column(5, #1=plain, 2=bold, 3=italic, 4=bold italic
                                                     selectInput(inputId = "numFont", label = "Font Style", 
                                                                 choices = c(Normal = 1, Bold = 2, Italic = 3, `Bold Italic` = 4), 
                                                                 selected = 2)
                                              )
                                            ),
                                            fixedRow(
                                              column(4, 
                                                     numericInput(inputId = "numOffY", label = "Vertical Offset", 
                                                                  value = 0.3, min = -1, max = 1, step = 0.01)
                                              ),
                                              column(4, 
                                                     numericInput(inputId = "numOffX", label = "Horizontal Offset", 
                                                                  value = 0, min = -1, max = 1, step = 0.01)
                                              ),
                                              column(4, 
                                                     numericInput(inputId = "numOffPos", label = "Position Offset", 
                                                                  value = 0, min = 0, max = 10000, step = 1)
                                              )
                                            ),
                                            tags$hr(),
                                            helpText("Font Color", style = subtitlesStyle),
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "numCol1", label = "Polar / Basic", value = "black")
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "numCol2", label = "Polar / Acidic", value = "black")
                                              )),
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "numCol3", label = "Polar / Uncharged", value = "black")
                                              ),
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "numCol4", label = "Nonpolar", value = "black")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     colourpicker:: colourInput(inputId = "numCol5", label = "Unknown Residue", value = "black")
                                              )
                                            )
                                   )
                                 )
                        ),
                        tabPanel("Guides",
                                 helpText("Wheel Guides", style = subtitlesStyle),
                                 fixedRow(
                                   column(4,
                                          selectInput(inputId = "showWheelGuide", label = "Show",
                                                      choices = c("Yes", "No"), selected = "No")
                                   ),
                                   column(4,
                                          selectInput(inputId = "wheelGuideLty", label = "Line Type",
                                                      choices = ltys, selected = "33")
                                   ),
                                   column(4,
                                          numericInput(inputId = "wheelGuideLwd", label = "Line Width",
                                                       value = 1, min = 0, max = 5, step = 1)
                                   )
                                 ),
                                 fixedRow(
                                   column(6,
                                          numericInput(inputId = "wheelGuideCol", label = "Gray Level",
                                                       value = 0.5, min = 0, max = 1, step = 0.05)
                                   ),
                                   column(6,
                                          selectInput(inputId = "showBox", label = "Show Outer Box",
                                                      choices = c("Yes", "No"), selected = "Yes")
                                   )
                                 ),
                                 tags$hr(),
                                 helpText("Net Guides", style = subtitlesStyle),
                                 fixedRow(
                                   column(4,
                                          selectInput(inputId = "showNetGuide", label = "Show",
                                                      choices = c("Yes", "No"), selected = "No")
                                   ),
                                   column(4,
                                          selectInput(inputId = "netGuideLty", label = "Line Type",
                                                      choices = ltys, selected = "33")
                                   ),
                                   column(4,
                                          numericInput(inputId = "netGuideLwd", label = "Line Width",
                                                       value = 1, min = 0, max = 5, step = 1)
                                   )
                                 ),
                                 fixedRow(
                                   column(6,
                                          numericInput(inputId = "netGuideCol", label = "Gray Level",
                                                       value = 0.5, min = 0, max = 1, step = 0.05)
                                   ),
                                   column(6,
                                          selectInput(inputId = "showBoxNet", label = "Show Outer Box",
                                                      choices = c("Yes", "No"), selected = "No")
                                   )
                                 )
                        ),
                        tabPanel("Legend",
                                 tabsetPanel(
                                   tabPanel("Position and Labels",
                                            helpText("Wheel", style = subtitlesStyle),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "showLeg", label = "Show",
                                                                 choices = c("Yes", "No"), selected = "Yes")
                                              ),
                                              column(4,
                                                     numericInput(inputId = "legX", label = "Relative x",
                                                                  value = 0.8, min = -2, max = 2, step = 0.01)
                                              ),
                                              column(4,
                                                     numericInput(inputId = "legY", label = "Relative y",
                                                                  value = 0.9, min = -2, max = 2, step = 0.01)
                                              )
                                            ),
                                            helpText("Net", style = subtitlesStyle),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "showLegNet", label = "Show",
                                                                 choices = c("Yes", "No"), selected = "Yes")
                                              ),
                                              column(4,
                                                     numericInput(inputId = "legXNet", label = "Relative x",
                                                                  value = 0, min = -10, max = 10, step = 0.05)
                                              ),
                                              column(4,
                                                     numericInput(inputId = "legYNet", label = "Relative y",
                                                                  value = -1, min = -15, max = 3, step = 0.05)
                                              )
                                            ),
                                            tags$hr(),
                                            helpText("Labels", style = subtitlesStyle),
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "leg1", label = "Polar / Basic",
                                                               value = "Polar / Basic")
                                              ),
                                              column(6,
                                                     textInput(inputId = "leg2", label = "Polar / Acidic",
                                                               value = "Polar / Acid")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,textInput(inputId = "leg3", label = "Polar / Uncharged",
                                                                 value = "Polar / Uncharged")
                                                     
                                              ),
                                              column(6,
                                                     textInput(inputId = "leg4", label = "Nonpolar",
                                                               value = "Nonpolar")
                                              )
                                            ),
                                            textInput(inputId = "leg5", label = "Group 5 Legend",
                                                      value = "Unknown Residue")
                                   ),
                                   tabPanel("Aesthetics",
                                            numericInput(inputId = "legCex", label = "Font Size",
                                                         value = 16, min = 1, max = 30, step = 1),
                                            tags$hr(),
                                            helpText("Pattern Density (Lines)", style = subtitlesStyle),
                                            helpText("Note that these may not display correctly on the online version because of the low DPI (72). Make sure to check the exported output."),
                                            fixedRow(
                                              column(4,
                                                     numericInput(inputId = "legDen1", "Polar / Basic",
                                                                  value = 20, min = 1, max = 50)
                                              ),
                                              column(4,
                                                     numericInput(inputId = "legDen2", "Polar / Acidic",
                                                                  value = 24, min = 1, max = 50)
                                              ),
                                              column(4,
                                                     numericInput(inputId = "legDen3", "Polar / Uncharged",
                                                                  value = 31, min = 1, max = 50)
                                              )
                                            ),
                                            fixedRow(
                                              column(4,
                                                     numericInput(inputId = "legDen4", "Nonpolar ",
                                                                  value = 20, min = 1, max = 50)
                                              ),
                                              column(4,
                                                     numericInput(inputId = "legDen5", "Group 5",
                                                                  value = 20, min = 1, max = 50)
                                              )
                                            )
                                   )
                                 )
                        ),
                        tabPanel("Groups",
                                 tabsetPanel(
                                   tabPanel("Residues",
                                            helpText("Use one-letter code"),
                                            helpText("Residues main group", style = subtitlesStyle),
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "grp1", label = "Polar / Basic",
                                                               value = "RHK")
                                              ),
                                              column(6,
                                                     textInput(inputId = "grp2", label = "Polar / Acidic",
                                                               value = "DE")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "grp3", label = "Polar / Uncharged",
                                                               value = "STNQC")
                                              ),
                                              column(6,
                                                     textInput(inputId = "grp4", label = "Nonpolar ",
                                                               value = "AGVILMFYWP")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "grp5", label = "Group 5",
                                                               value = "X", width = "165px")
                                              ),
                                              column(6,
                                                     actionButton("grpReset", "Reset to Default", style = "margin-top: 25px")
                                              )
                                            ),
                                            helpText("Residues interactions (Net)", style = subtitlesStyle),
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "grpNonpolar", label = "Nonpolar",
                                                               value = "VILMFYW")
                                              ),
                                              column(6,
                                                     textInput(inputId = "grpHydro", label = "Hydrogen Bond",
                                                               value = "STNQY")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "grpAcid", label = "Acid",
                                                               value = "DE")
                                              ),
                                              column(6,
                                                     textInput(inputId = "grpBasic", label = "Basic",
                                                               value = "RHK")
                                              )
                                            ),
                                            actionButton("grpBondReset", "Reset to Default")
                                   ),
                                   tabPanel("Labels",
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "grp1Lab", label = "Group 1",
                                                               value = "Polar / Basic")
                                              ),
                                              column(6,
                                                     textInput(inputId = "grp2Lab", label = "Group 2",
                                                               value = "Polar / Acidic")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "grp3Lab", label = "Group 3",
                                                               value = "Polar / Uncharged")
                                              ),
                                              column(6,
                                                     textInput(inputId = "grp4Lab", label = "Group 4",
                                                               value = "Nonpolar")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     textInput(inputId = "grp5Lab", label = "Group 5",
                                                               value = "Unkown Residue", width = "165px")
                                              ),
                                              column(6,
                                                     actionButton("grpUpdUi", "Update Interface", style = "margin-top: 25px")
                                              )
                                            ),
                                            fixedRow(
                                              column(6,
                                                     actionButton("grpLabToLeg", "Copy to Legend", style = "margin-top: 25px")
                                              ),
                                              column(6,
                                                     actionButton("grpLabDefault", "Reset to defaults", style = "margin-top: 25px")
                                              )
                                            )
                                   )
                                 )
                                 
                        ),
                        tabPanel("Title",
                                 #helpText("Adjust the xy postions and top-margin using the 'Wheel' or 'Net' tabs if necessary"),
                                 fixedRow(
                                   column(4,
                                          selectInput(inputId = "showTitle", label = "Show",
                                                      choices = c("Yes", "No"), selected = "No")
                                   ),
                                   column(4,
                                          actionButton(inputId = "titleSeq", label = "Copy from Sequence",
                                                       style = "margin-top: 25px")
                                   )
                                 ),
                                 fixedRow(
                                   column(12,
                                          textInput(inputId = "txTitle", label = "Text",
                                                    value = "Title", width = "100%")
                                   )
                                 ),
                                 hr(),
                                 helpText("Wheel Positioning"),
                                 fixedRow(
                                   column(4,
                                          numericInput(inputId = "xTitle", label = "Relative x", 
                                                       value = 0, min = -10, max = 20, step = 0.1)
                                   ),
                                   column(4,
                                          numericInput(inputId = "yTitle", label = "Relative y", 
                                                       value = 0, min = -10, max = 20, step = 0.1)
                                   ),
                                   column(4,
                                          numericInput(inputId = "cexTitle", label = "Font Size", 
                                                       value = 20, min = 1, max = 30, step = 1)
                                   )
                                 ),
                                 hr(),
                                 helpText("Net Positioning"),
                                 fixedRow(
                                   column(4,
                                          numericInput(inputId = "xTitleNet", label = "Relative x", 
                                                       value = 2.3, min = -100, max = 300, step = 0.1)
                                   ),
                                   column(4,
                                          numericInput(inputId = "yTitleNet", label = "Relative y", 
                                                       value = 10.6, min = -100, max = Inf, step = 0.1)
                                   ),
                                   column(4,
                                          numericInput(inputId = "cexTitleNet", label = "Font Size", 
                                                       value = 14, min = 1, max = 30, step = 1)
                                   )
                                 ),
                                 actionButton("netAutoMar", "Auto Adjust Position")
                        ),
                        tabPanel("Net",
                                 tabsetPanel(
                                   tabPanel("General",
                                            fixedRow(
                                              column(4,
                                                     actionButton("alphaHelix", "alpha helix")
                                              ),
                                              column(4,
                                                     actionButton("three10Helix", "3_10 helix")
                                              ),
                                              column(4,
                                                     actionButton("piHelix", "pi helix")
                                              )
                                            ),
                                            fixedRow(
                                              column(4,
                                                     numericInput("netPerTurn", "Res. per turn",
                                                                  value=3.6, min=2, max=8.0, step=0.01)
                                              ),
                                              column(4,
                                                     selectInput(inputId = "netDirection", label = "Direction",
                                                                 choices = c(`Right-Left` = 1, `Left-Right` = -1),
                                                                 selected = 1)
                                              ),
                                              column(4,
                                                     numericInput("netStartOff", "Start Offset",
                                                                  value=0, min=0, max=10, step=0.05)
                                              )
                                            ),
                                            #                                             fixedRow(
                                            #                                               column(4,
                                            #                                                      numericInput("netDiameter", "Diameter (Å)",
                                            #                                                                   value=4.6, min=3.0, max=6.0, step=0.1)
                                            #                                                      ),
                                            #                                               column(3,
                                            #                                                      numericInput("netPitch", "Pitch (Å)",
                                            #                                                                   value=5.4, min=4.5, max=6.5, step=0.01)
                                            #                                               ),
                                            #                                               column(5,
                                            #                                                      numericInput("netTrans", "Translations (Å)",
                                            #                                                                   value=1.5, min=1.0, max=2.5, step=0.01)
                                            #                                               )
                                            #                                             ),
                                            fixedRow(
                                              column(6,
                                                     numericInput("netPadTop", "Padding Top",
                                                                  value = 0, min = 0, max = 3, step = 0.05)
                                              ),
                                              column(6,
                                                     numericInput("netPadBot", "Padding Bottom",
                                                                  value = 0.2, min = 0, max = 1, step = 0.01)
                                              )
                                            ),
                                            fixedRow(
                                              column(4,
                                                     numericInput("netPadL", "Padding Left",
                                                                  value = 0, min = 0, max = 1, step = 0.01)
                                              ),
                                              column(4,
                                                     numericInput("netPadR", "Padding Right",
                                                                  value = 0, min = 0, max = 1, step = 0.01)
                                              ),
                                              column(4,
                                                     selectInput(inputId = "netShowLimits", label = "Show Border",
                                                                 choices = c("Yes", "No"), selected = "Yes")
                                              )
                                            ),
                                            helpText("Change this to make the polygons proportional and adjust the figre size, if necessary:"), 
                                            fixedRow(column(6,
                                                            numericInput("netProp", "Proportion Factor",
                                                                         value = 2.8, min = 0.05, max = 10, step = 0.05)
                                            ),
                                            column(6,
                                                   numericInput("netWidth", "Figure Extra Width",
                                                                value = 240, min = 0, max = 500, step = 5)
                                            ))
                                            
                                   ),
                                   tabPanel("Interactions",
                                            selectInput(inputId = "netShowInteractions", label = "Show Interactions",
                                                        choices = c("Yes", "No"), selected = "Yes"),
                                            helpText("Nonpolar interactions"),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "bond1Ty", label = "Type",
                                                                 choices = ltys, selected = 1)
                                              ),
                                              column(3,
                                                     numericInput(inputId = "bond1Wd", label = "Width",
                                                                  value = 3, min = 0, max = 10, step = 1)
                                              ),
                                              column(5,
                                                     colourpicker:: colourInput(inputId = "bond1Col", label = "Color",
                                                                          value = "black")
                                              )
                                            ),
                                            helpText("Acid-base interactions"),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "bond2Ty", label = "Type",
                                                                 choices = ltys, selected = "33")
                                              ),
                                              column(3,
                                                     numericInput(inputId = "bond2Wd", label = "Width",
                                                                  value = 3, min = 0, max = 10, step = 1)
                                              ),
                                              column(5,
                                                     colourpicker:: colourInput(inputId = "bond2Col", label = "Color",
                                                                          value = "black")
                                              )
                                            ),
                                            helpText("Hydrogen bonds"),
                                            fixedRow(
                                              column(4,
                                                     selectInput(inputId = "bond3Ty", label = "Type",
                                                                 choices = ltys, selected = "12")
                                              ),
                                              column(3,
                                                     numericInput(inputId = "bond3Wd", label = "Width",
                                                                  value = 3, min = 0, max = 10, step = 1)
                                              ),
                                              column(5,
                                                     colourpicker:: colourInput(inputId = "bond3Col", label = "Color",
                                                                          value = "black")
                                              )
                                            )
                                   )
                                 )
                        ),
                        tabPanel("Wheel",
                                 helpText("Connection Lines", style = subtitlesStyle),
                                 fixedRow(
                                   column(4,
                                          numericInput(inputId = "conLineWd", label = "Width", 
                                                       value = 6, min = 0, max = 20, step = 1)
                                   ),
                                   column(4, 
                                          numericInput(inputId = "maxlinegray", label = "Max Gray", 
                                                       value = 0.7, min = 0, max = 1, step = 0.01)
                                   ),
                                   column(4, 
                                          numericInput(inputId = "minlinegray", label = "Min Gray", 
                                                       value = 0.1, min = 0, max = 1, step = 0.01)
                                   )
                                 ),
                                 tags$hr(),
                                 helpText("Wheel Properties", style = subtitlesStyle),
                                 fixedRow(
                                   column(6,
                                          numericInput(inputId = "innRatio", label = "Spiral Factor", 
                                                       value = 0, min = 0, max = 1, step = 0.01)
                                   ),
                                   column(6,
                                          numericInput(inputId = "nestRatio", label = "Nesting Factor", 
                                                       value = 0.5, min = 0, max = 2, step = 0.01)
                                   )
                                 ),
                                 fixedRow(
                                   column(6,
                                          selectInput(inputId = "direction", label = "Direction", 
                                                      choices = c("Clockwise", "Counterclockwise"),
                                                      selected = "Clockwise")
                                   ),
                                   column(3,
                                          numericInput(inputId = "period", label = "Period", 
                                                       value = 18, min = 5, max = 30, step = 1)
                                   ),
                                   column(3,
                                          numericInput(inputId = "perStep", label = "Step", 
                                                       value = 4, min = 1, max = 10, step = 1)
                                   )
                                 ),
                                 helpText("Figure Size", style = subtitlesStyle),
                                 fixedRow(
                                   column(4,
                                          numericInput(inputId = "figWd", label = "Figure Width",
                                                       value = 775, min = 400, max = 1000, step = 5)
                                   ),
                                   column(4,
                                          numericInput(inputId = "marR", label = "Padding-Right",
                                                       value = 12, min = 0, max = 30, step = 0.1)
                                   ),
                                   column(4,
                                          numericInput(inputId = "marT", label = "Padding-Top",
                                                       value = 0, min = 0, max = 10, step = 0.1)
                                   )
                                 ),
                                 helpText("Use Width = 600 and Padding = 0 for a square image without legend.")
                        ),
                        tabPanel("Styles",
                                 helpText("Some pre-defined styling actions"),
                                 helpText("Residues Representation", style = subtitlesStyle),
                                 helpText("Colors"),
                                 fixedRow(
                                   column(3,
                                          actionButton("reset_input2", "Default")
                                   ),
                                   column(3,
                                          actionButton("style_terrain", "Terrain")
                                   ),
                                   column(3,
                                          actionButton("style_bw", "B & W")
                                   ),
                                   column(3,
                                          actionButton("style_gray", "Grayscale")
                                   )
                                 ),
                                 helpText("Shapes"),
                                 fixedRow(
                                   column(3,
                                          actionButton("style_circle", "Circles")
                                   ),
                                   column(3,
                                          actionButton("style_squares", "Squares")
                                   ),
                                   column(2,
                                          actionButton("style_shapes", "All")
                                   ),
                                   column(4,
                                          actionButton("style_pattern", "Circle Pattern")
                                   )
                                 ),
                                 actionButton("style_labdown", "Switch Label Style"),
                                 helpText("Wheel", style = subtitlesStyle),
                                 fixedRow(
                                   column(6,
                                          actionButton("style_legend", "Toggle Wheel Legend")
                                   ),
                                   column(6,
                                          actionButton("style_round", "Toggle Spiral")
                                   )
                                 ),
                                 helpText("Net", style = subtitlesStyle),
                                 fixedRow(
                                   column(4,
                                          actionButton("style_netLeg0", "No Legend")
                                   ),
                                   column(4,
                                          actionButton("style_netLeg1", "Legend Above")
                                   ),
                                   column(4,
                                          actionButton("style_netLeg2", "Legend Beside")
                                   )
                                 )
                        ),
                        tabPanel("Tools",
                                 tags$br(),
                                 downloadButton("bt_export", "Export Settings"),
                                 #actionButton("bt_rdata", "Save RData"),
                                 #actionButton("bt_rdataload", "Load RData"),
                                 tags$br(),
                                 tags$br(),
                                 fileInput("fileImport", "Import Style File",
                                           accept=c(".hws")),
                                 tags$script('
        Shiny.addCustomMessageHandler("resetFileInputHandler", function(x) {   
          var el = $("#" + x);
          el.replaceWith(el = el.clone(true));
          //var id = "#" + x + "_progress";     
          //$(id).css("visibility", "hidden");
        });
      ')
                                 #actionButton("do_import", "Activate Import")
                        ),
                        tabPanel("Export",
                                 helpText("Export Options", style = subtitlesStyle),
                                 selectInput(inputId = "expFormat", label = "Format",
                                             choices = c("PNG", "TIFF", "PDF", "JPG"), selected = "PNG"),
                                 selectInput(inputId = "expDPI", label = "DPI",
                                             choices = c(72, 150, 300, 600), selected = "72"),
                                 fixedRow(
                                   column(6,
                                          downloadButton("plotDown", "Export Wheel")),
                                   column(6,
                                          downloadButton("plotDownNet", "Export Net")
                                   )
                                 )
                        )
                      )
             )
      ),
      column(8,
             tabsetPanel(
               tabPanel("Wheel",
                        checkboxInput("autoWheel", "Automatic Preview", value = FALSE),
                        #conditionalPanel(condition = "input.autoWheel == true",
                        plotOutput(outputId = "helicalPlot", height = "600px", width = "600px")#,
                        #                         tags$hr(),
                        #                         helpText("Export Options", style = subtitlesStyle),
                        #                         fixedRow(
                        #                           column(4,
                        #                                  selectInput(inputId = "expFormat", label = "Format",
                        #                                              choices = c("PNG", "TIFF (lzw)", "PDF", "JPG"), selected = "PNG")
                        #                           ),
                        #                           column(4,
                        #                                  selectInput(inputId = "expDPI", label = "DPI",
                        #                                              choices = c(72, 150, 300, 600), selected = "72")
                        #                           ),
                        #                           column(4,
                        #                                  downloadButton("plotDown", "Save Image")
                        #                           ),
                        #                           tags$style(type='text/css', "#plotDown { margin-top: 25px;}")
                        #                         )
                        #)
               ),
               tabPanel("Net",
                        checkboxInput("autoNet", "Automatic Preview", value = FALSE),
                        plotOutput(outputId = "netPlot", 
                                   height = "600px", width = "300px")
               ),
               tabPanel("About",
                        tags$br(),
                        tags$p("Developed by Alan R. Mól, Wagner Fontes, Mariana S. Castro."),
                        HTML("<p><a href=\"http://lbqp.unb.br\" target=\"_blank\">Laboratório de Bioquímica e Química de Proteínas (LBQP)</a> - Universidade de Brasília - Brazil</p>"),
                        tags$p("This application is available free of charge. If you would like to cite it, please use the following reference:"),
                        HTML("<p><a href=\"https://doi.org/10.26502/jbsb.5107082\" target=\"_blank\">Mól AR, Castro MS, Fontes W. NetWheels: A Web Application to Create High Quality Peptide Helical Wheel and Net Projections. Journal of Bioinformatics and Systems Biology. 7 (2024): 98-100; doi: https://doi.org/10.1101/416347</a></p>"),
                        HTML("<p>Source code, local installation instructions and other information can be found on <a href=\"https://github.com/molx/NetWheels\">Github</a>."),
                        tags$br()#,
                        #tags$p("To-do list:"),
                        #tags$ul(tags$li("Enable polygon size based on physical properties of the residues, like originally used by Dunnill (1968)"),
                               # tags$li("Add legend-support for different shapes. Currently, only colors and patterns are considered.")
                        #)
               ),
               tabPanel("Help",
                        tags$br(),
                        tags$p("For quick examples of styles and also one-click changes, check the 'Styles' tab. Advanced ajustments, as well as other features, are available using the other tabs, explained below:"),
                        # The tabs are created in the do.call call, to avoid repeating the expression
                        do.call(tabsetPanel, lapply(seq_along(tabs), function(hT) {
                          tabPanel(tabs[hT],
                                   tags$br(),
                                   # Each line from the .txt file is transformed into a p tag. 
                                   lapply(seq_along(helpTexts[[hT]]), function(h) tags$p(helpTexts[[hT]][h], style="text-indent: 50px;"))
                          )
                        }))
               )
             )
      )
    ),
    tags$br(),
    tags$br(),
    tags$br(),
    tags$br(),
    tags$div(
      tags$footer("Developed by Alan R. Mól, Wagner Fontes, Mariana S. Castro, Universidade de Brasília - Brazil")
    )
  )
  
)
)
