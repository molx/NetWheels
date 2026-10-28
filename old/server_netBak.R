library(shiny)

shinyServer(function(input, output, session) {
  
  # All aminoacids and their classification
  amin <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "G", "P", "A", "V", "I", "L", "M", "F", "Y", "W", "X")
  hb   <- c( 1 ,  1 ,  1 ,  2 ,  2 ,  3 ,  3 ,  3 ,  3 ,  3 ,  3 ,  3 ,  3 ,  4 ,  4,   4 ,  4 ,  4 ,  4 ,  4 ,  4  , 5 )
  #        charged/basic    -acid-    ------------polar/uncharged---------    -------------hydrophibic-------------
  
  # The 'X' aminoacid represents any aminoacid. All characters not present in amin will be replaced with X, unless a 3-letter format is detected
  
  amin3 <-c("Arg", "His", "Lys", "Asp", "Glu", "Ser", "Thr", "Asn", "Gln", "Cys", "Sec", "Gly", "Pro", "Ala",
            "Val", "Ile", "Leu", "Met", "Phe", "Tyr", "Trp")
  
  # Peptide groups in terms of bonds, some do miss
  
  apolar <- c("V", "I", "L", "M", "F", "Y", "W")
  acid <- c("D", "E")
  basic <- c("R", "H", "K")
  hydrobond <- c("S", "N", "Q")
  
  addPolygon <- function(shape, x, y, size, col, lwd, border) {
    if (shape == "Triangle") {
      a <- 2*size*2/sqrt(3) # Length of triangle size. r.res is half of the height.
      xi <- c(x, x + a/2, x + a) - a/2 
      yi <- c(y, y + size*2, y) - size*0.8 # Adjusts ypos to avoid overlap of label and triangle thin top
    } else if (shape == "Square") {
      xi <- c(x, x, x + size*2, x + size*2) - size
      yi <- c(y, y + size*2, y + size*2, y) - size
    } else if (shape == "Diamond") {
      d <- 0.85 * size*2*sqrt(2) # Diagonal of the square, since it's rotated. 0.8 lowers it because the vertices make it look bigger
      xi <- c(x, x + d/2, x + d, x + d/2) - d/2
      yi <- c(y, y + d/2, y, y - d/2)
    } else if (shape == "Hexagon") {
      d <- 2*size*2/sqrt(3) # Hexagon diagonal
      l <- 2*size/sqrt(3) # Distance between two vertices
      p1 <- (d - l)/2 # y-axis movement step
      xi <- c(x, x - size, x - size, x, x + size, x + size)
      yi <- c(y, y + p1, y + p1 + l, y + d, y + p1 + l, y + p1) - d/2
    } else { #Circle
      thetas <- seq(0, 2*pi, length.out = 100)
      xi <- size * sin(thetas) + x
      yi <- size * cos(thetas) + y
    }
    polygon(xi, yi, col = col, lwd = lwd,
            border = border)
  }
  
  rotvec <- function (v, n) {
    # Taken from wavethresh::guyrot
    l <- length(v)
    n <- n%%l
    if (n == 0) 
      return(v)
    tmp <- v[(l - n + 1):l]
    v[(n + 1):l] <- v[1:(l - n)]
    v[1:n] <- tmp
    v
  }
  
  #### Code for Helical Net ####

  getHeight <- function() {
    pepSeq <- input$seq
    n <- nchar(pepSeq)
    rws <- ceiling(n/7)
    height <- 200*rws*input$ydist
    height
  }

  netPlot <- function() {
    pepSeq <- input$seq
    
    n <- nchar(pepSeq)
    
    width <- 300
    rws <- ceiling(n/7)
    #ydist <- 1
    height <- 200*rws*input$ydist
    
    imarx <- width*input$imarx # margem interna x
    imary <- height*input$imary # margem interna y
    par(mar=c(0,0,0,0) + 0.1)
    plot(1, type = "n", xlim = c(0, width), ylim = c(0, height),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         xaxs = "i", yaxs = "i", frame.plot = input$showBoxNet == "Yes")
    #abline(h=c(imary, height-imary), v = c(imarx, width-imarx), lty = 2)
    
    lfrac <- 0.15
    l <- width*lfrac
    if (input$netDirection == "1") {
      xpos <- c(4:1, 3:1+0.5)*(width - 2*imarx)/5 + imarx 
    } else {
      xpos <- c(1:4, 1:3+0.5)*(width - 2*imarx)/5 + imarx 
    }
    
    # Changing start positions:
    xpos <- rotvec(xpos, as.numeric(input$netStart) - 1)
    
    ptx <- rep(xpos, ceiling(n/7))[1:n]
    
    pty <- rev((1:n)*(height - 2*imary)/(n+1)) + imary
    
    # Adding guides behind
    
    if (input$showNetGuide == "Yes") {
      ltyNetGuide <- as.numeric(input$netGuideLty)
      lwdNetGuide <- as.numeric(input$netGuideLwd)
      colNetGuide <- gray(1-as.numeric(input$netGuideCol))
      
      print(c(ltyNetGuide, lwdNetGuide, colNetGuide))
      
      lines(ptx, pty, lty = ltyNetGuide, lwd = lwdNetGuide, col = colNetGuide) 
    }
    
    res <- unlist(strsplit(pepSeq, ""))
    fills <- hb[match(res, amin)]
    
    shp <- c(input$shp1, input$shp2, input$shp3, input$shp4, input$shp5)
    cores <- c(input$col1, input$col2, input$col3, input$col4, input$col5)
    bordersShow <- c(input$circBorder1, input$circBorder1,
                     input$circBorder3, input$circBorder4, input$circBorder5)
    
    if (any(bordersShow == "Yes")) {
      bordersCol <- c(input$circBorderCol1, input$circBorderCol2,
                      input$circBorderCol3, input$circBorderCol4, input$circBorderCol5)
      bordersCol[bordersShow != "Yes"] <- NA
      bordersCol <- bordersCol[fills]
      bordersWd <- as.numeric(c(input$circBorderWd1, input$circBorderWd2,
                                input$circBorderWd3, input$circBorderWd4, input$circBorderWd5))
      bordersWd <- bordersWd[fills]
    } else {
      bordersCol <- rep(NA, length(nMin))
      bordersWd <- rep(1, length(nMin))
    }
    
    shapes <- shp[fills]
    
    nres <- nchar(pepSeq) 
    numLabs <- c(input$labCol1, input$labCol2, input$labCol3, input$labCol4, input$labCol5)
    
    # Making labels
    labType <- as.numeric(input$labType)
    
    resLab <-  if (labType == 0) {
      ""
    } else if (labType == 1) {
      res
    } else if (labType == 2) {
      amin3[match(res, amin)]
    } else if (labType == 3) {
      paste0(res, seq_len(nres))
    } else if (labType == 4) {
      seq_len(nres)
    }
    
    bondWd <- c(apolar = input$bond1Wd, acba = input$bond2Wd, hydro = input$bond3Wd)
    bondTy <- c(apolar = input$bond1Ty, acba = input$bond2Ty, hydro = input$bond3Ty)
    bondCol <- c(apolar = input$bond1Col, acba = input$bond2Col, hydro = input$bond3Col)
    
    for (i in seq_len(n)) {
      # Detecting Bonds 
      
      bond3 <- TRUE
      if (res[i] %in% apolar && res[i + 3] %in% apolar) {
        bond3type <- "apolar"
      } else if ((res[i] %in% basic && res[i + 3] %in% acid) ||
                 (res[i] %in% acid && res[i + 3] %in% basic)) {
        bond3type <- "acba"
      } else if (res[i] %in% hydrobond && res[i + 3] %in% hydrobond) {
        bond3type <- "hydro"
      } else {
        bond3 <- FALSE
        bond3type <- "none"
      } 
      
      bond4 <- TRUE
      if (res[i] %in% apolar && res[i + 4] %in% apolar) {
        bond4type <- "apolar"
      } else if ((res[i] %in% basic && res[i + 4] %in% acid) ||
                 (res[i] %in% acid && res[i + 4] %in% basic)) {
        bond4type <- "acba"
      } else if (res[i] %in% hydrobond && res[i + 4] %in% hydrobond) {
        bond4type <- "hydro"
      } else {
        bond4 <- FALSE
        bond4type <- "none"
      } 
      # Adding bonds
      #if (i < n-3 && res[i] == "I" && res[i+3] == "I") {
      
      if (input$netDirection == "1") {
        # The min/max choice for broken bonds depends on the direction of the net
        # We select the min or max functions to be the extremes depending on the direction
        ext1 <- max
        ext2 <- min
      } else {
        ext1 <- min
        ext2 <- max
      }
      
      lty3 <- if(nchar(bondTy[bond3type]) > 1) bondTy[bond3type] else as.numeric(bondTy[bond3type])
      lty4 <- if(nchar(bondTy[bond4type]) > 1) bondTy[bond4type] else as.numeric(bondTy[bond4type])
      
      if (bond3) {
        if(ptx[i] == ext1(ptx)) {
          xdif <- (ptx[i] - ptx[i+3])/5
          xp2 <- ptx[i] + xdif/2
          yp2 <- pty[i] - (pty[i] - pty[i+3]+l)/3
          segments(ptx[i+3], pty[i+3]+l/2, ext2(ptx) - xdif/2, yp2,
                   lwd = as.numeric(bondWd[bond3type]), lty = lty3,
                   col = bondCol[bond3type])
        } else {
          xp2 <- ptx[i+3]
          yp2 <- pty[i+3]+l/2
        }
        segments(ptx[i], pty[i]-l/2, xp2, yp2,
                 lwd = as.numeric(bondWd[bond3type]), lty = lty3,
                 col = bondCol[bond3type])
      }
      #if (i < n-4 && res[i] == "I" && res[i+4] == "I") {
      if (bond4) {
        if(ptx[i] == ext2(ptx)) {
          xdif <- (ptx[i+4] - ptx[i])/5
          xp2 <- ptx[i] - xdif
          yp2 <- pty[i] - (pty[i] - pty[i+4]+l/2)/2
          segments(ptx[i+4], pty[i+4]+l/2, ext1(ptx) + xdif/2, yp2,
                   lwd = as.numeric(bondWd[bond4type]), lty = lty4,
                   col = bondCol[bond4type])
        } else {
          xp2 <- ptx[i+4]
          yp2 <- pty[i+4]+l/2
        }
        segments(ptx[i], pty[i]-l/2, xp2, yp2,
                 lwd = as.numeric(bondWd[bond4type]), lty = lty4,
                 col = bondCol[bond4type])
      }
      
      # Adding Residues Polygons
      bwd <- bordersWd[i]
      addPolygon(shape = shapes[i], 
                 x = ptx[i], y = pty[i],
                 size = l/2, col = cores[fills][i],
                 lwd = bwd, border = bordersCol[i])
#       text(ptx[i], pty[i], label = resLab[i],
#            cex = input$labCex/10, font = as.numeric(input$labFont),
#            col = numLabs[fills])
    }
    
    text(ptx, pty, label = resLab,
         cex = input$labCex/10, font = as.numeric(input$labFont),
         col = numLabs[fills])
    
    # Adding title/name. Later to be on top
    
    if (input$showTitle == "Yes") {
      main <- input$txTitle
      cex.main <- input$cexTitleNet/10
      text(x = input$xTitleNet, y = input$yTitleNet, label = main,
           cex = cex.main, adj = c(0.5, 0.5), font = 2)
    }

  }
  
  output$netPlot <- renderPlot(
    # The auto option is disbled by default to avoid unecessary plot creation 
    # when another projection is being tested, since this is rather slow
    # and require considerable resources when online
    {
      if(input$autoNet || is.null(lastWheel)) {
        netPlot()
        lastNet <<- recordPlot()
      } else {
        replayPlot(lastNet)
      }
    },
    height = getHeight, width = 300)
  
  
  #### Code for Helical Wheel ####
  
  # Function to add stripes inside circles
  
  circStripes <- function(r, n, x=0, y=0, col = "black", angle = 0) {
    h <- seq(0, r, length.out = n - ((n-1) %/% 3))
    angs <- 2 * acos((r-h)/r)
    chord <- r * sqrt(2 - 2*cos(angs))
    hp <- c(h-r, r-h)
    x1 <- c(-chord/2+x, -chord/2+x, +chord/2+x, +chord/2+x)
    y1 <- c(hp + y, hp + y)
    a <- angle*pi/180
    
    x2 <- (x1 - x) * cos(a) - (y1 - y) * sin(a) + x
    y2 <- (x1 - x) * sin(a) + (y1 - y) * cos(a) + y
    len <- length(x1)
    segments(x2[1:(len/2)], y2[1:(len/2)], x2[(len/2+1):len], y2[(len/2+1):len], col = col)
  }
  
  # Wheel radius. Changing this just messes with some relative positions offsets.
  r <- 10 #input$wheelsize # Raio do circulo grande  da figura
  
  lastWheel <- NULL
  
  helicalPlot <- function() {
    
    nHel <- input$period
    perStep <- input$perStep
    
    #pos <- c(1, 6, 11, 16, 3, 8, 13, 18, 5, 10, 15, 2, 7, 12, 17, 4, 9, 14) 
    
    # Creating the minimum positions for 18 residues
    posMin <- 1
    
    for (i in 2:nHel) {
      y <- posMin[i-1] + perStep + 1
      if (y > nHel) y <- y - nHel
      posMin[i] <- y
    }
    
    # Reading the pepseq
    pepSeq <- input$seq 
    
    if ((nchar(gsub("[A-Z]", "", pepSeq))/nchar(gsub("[a-z]", "", pepSeq))) == 2) {
      # Tries to detect 3-letter code and convert to 1-letter
      pepSeq <- amin[match(pepSeq, amin3)]
    } else {
      # Converts everything to upper case. Making sure everything is matched and 
      # Also allows for unkown pep detection below safely
      pepSeq <- toupper(pepSeq)
      #updateTextInput(session, "seq", value = pepSeq)
    }
    
    # Replacing every character that isn't in 'amin' with 'X'
    splitSeq <- unlist(strsplit(pepSeq, ""))
    anyAmin <- which(!splitSeq %in% amin)
    if (length(anyAmin) > 0) {
      splitSeq[anyAmin] <- "X"
      pepSeq <- paste0(splitSeq, collapse = "")
      #updateTextInput(session, "seq", value = pepSeq)
    }
    
    # Number of aminoacids residues
    nres <- nchar(pepSeq) 
    
    output$resNumber <- renderUI(helpText(paste("Number of residues:", nres)))
    
    nHelicals <- ceiling(nres/nHel)
    
    angs <- seq(acos(0), acos(0)+2*pi, length.out = nHel + 1) #Angulos para posicao de cada circulo pequeno
    # acos(0) garante que o primeiro angulo terá como resultado um circulo no eixo x = 0 e y = max(y)
    # acos(0)+2*pi faz o circulo dar uma volta completa
    # nres + 1 é utilizado pois o ultimo ponto sempre coincide com o primeiro. Criamos um a mais e depois removemos o excedente.
    
    nMin <- if(nres < 18) 1:18 else 1:nres
    
    angs <- rep(angs, nHelicals)[nMin]
    
    # Order of the circles. 
    pos <- rep(posMin, nHelicals)[nMin]
    
    perim <- 2*pi*r #Perimetro do circulo grande
    
    innerOff <- rep(seq(1, 1-input$innRatio, length.out = nres)[order(pos[1:nHel])], nHelicals)[nMin]
    
    rotDirection <- if(input$direction == "Clockwise") -1 else 1
    
    ptx <- innerOff * (rotDirection*r) * cos(angs) # Pontos x
    
    pty <- innerOff * r * sin(angs) # Pontos y
    
    # Creating factor to multiply x and y and make them smaller for other circles
    
    nestFact <- (1/((nMin-1)%/%nHel + 1))^(1-input$nestRatio)
    
    ptxOrd <- ptx[pos] * nestFact
    
    ptyOrd <- pty[pos] * nestFact
    
    res.ratio <- input$circprop # Razao entre o tamanho (raio) do menor circulo em relacao ao maior (1)
    
    r.res.prop <- seq(res.ratio, 1, length.out = length(nMin)) # Measure that will be used to define the size of the residue polygon
    
    r.res.prop <- sqrt(r.res.prop / pi) # The circles should have the area, not the radius, proportional do the position. This fixes it based on A = pi * r * r
    
    # l.res.prop <- 
    
    res.space <- input$circsize # Tentativa: espaçamento entre circulos
    
    r.res <- rev((perim / sum(r.res.prop)) * r.res.prop * res.space) # Tentativa: raio dos circulos em função do raio total
    
    cores <- c(input$col1, input$col2, input$col3, input$col4, input$col5)
    
    res <- unlist(strsplit(pepSeq, ""))
    
    #par(mar = c(0,0,0,input$marR)+0.1, xpd = TRUE)
    #mpar()
    
    par(mar = c(0, 0, input$marT, input$marR)+0.1, xpd = TRUE)
    
    plot(ptx, pty, xlim = c(min(ptx, na.rm = TRUE) - max(r.res, na.rm = TRUE),
                            max(ptx, na.rm = TRUE) + max(r.res, na.rm = TRUE)),
         ylim = c(min(pty, na.rm = TRUE) - max(r.res, na.rm = TRUE),
                  max(pty, na.rm = TRUE) + max(r.res, na.rm = TRUE)),
         type = "n", xaxt = "n", yaxt = "n", frame.plot = FALSE,
         xlab = "", ylab = "", main = NA) # Preparando a area de plotagem, nao vai aparecer nada...
    
    
    # Adding circles guides. Before other stuff to stay behind everything
    
    if (input$showWheelGuide == "Yes") {
      angs.guide <- seq(0, 2*pi, length.out = nHel)
      
      ptx.guide <- - r * sin(angs.guide)
      pty.guide <- r * cos(angs.guide)
      
      lty <- as.numeric(input$wheelGuideLty)
      lwd <- as.numeric(input$wheelGuideLwd)
      col <- gray(1-as.numeric(input$wheelGuideCol))
      
      for (nest in unique(nestFact)) {
        lines(ptx.guide*nest, pty.guide*nest, lty = lty,
              lwd = lwd, col = col)  
      }
    }
    
    line_grad <- seq(1-input$maxlinegray, 1-input$minlinegray, length.out = nres)
    
    for (i in rev(seq_len(nres-1))) { #Plotando as conexões. 
      # De tras para frente para sobreposicao ficar correta
      sx0 <- ptxOrd[i]
      sx1 <- ptxOrd[i+1]
      sy0 <- ptyOrd[i]
      sy1 <- ptyOrd[i+1]
      segments(sx0, sy0, sx1, sy1, lwd = input$conLineWd, col =  gray(line_grad[i]))
    }
    
    # This matches which circle is of which category to properly select the colors and grid styles
    fills <- hb[match(res, amin)]
    
    bordersShow <- c(input$circBorder1, input$circBorder1,
                     input$circBorder3, input$circBorder4, input$circBorder5)
    
    if (any(bordersShow == "Yes")) {
      bordersCol <- c(input$circBorderCol1, input$circBorderCol2,
                      input$circBorderCol3, input$circBorderCol4, input$circBorderCol5)
      bordersCol[bordersShow != "Yes"] <- NA
      bordersCol <- bordersCol[fills]
      bordersWd <- as.numeric(c(input$circBorderWd1, input$circBorderWd2,
                                input$circBorderWd3, input$circBorderWd4, input$circBorderWd5))
      bordersWd <- bordersWd[fills]
    } else {
      bordersCol <- rep(NA, length(nMin))
      bordersWd <- rep(1, length(nMin))
    }
    
    # Ploting the polygons

    shp <- c(input$shp1, input$shp2, input$shp3, input$shp4, input$shp5)
    
    shapes <- shp[fills]
    
    for (i in nMin) {
      bwd <- bordersWd[i]
      if (is.na(bwd)) next # Avoiding errors with invalid lwd. NA throwns an error, skip because there's no circle to plot
      
      addPolygon(shape = shapes[i], 
                 x = ptxOrd[i], y = ptyOrd[i],
                 size = r.res[i], col = cores[fills][i],
                 lwd = bwd, border = bordersCol[i])
      
      # if (shapes[i] == "Triangle") { 
      #         a <- 2*r.res[i]*2/sqrt(3) # Length of triangle size. r.res is half of the height.
      #         xi <- c(ptxOrd[i], ptxOrd[i] + a/2, ptxOrd[i] + a) - a/2 
      #         yi <- c(ptyOrd[i], ptyOrd[i] + r.res[i]*2, ptyOrd[i]) - r.res[i]*0.8 # Adjusts ypos to avoid overlap of label and triangle thin top
      #       } else if (shapes[i] == "Square") { 
      #         xi <- c(ptxOrd[i], ptxOrd[i], ptxOrd[i] + r.res[i]*2, ptxOrd[i] + r.res[i]*2) - r.res[i]
      #         yi <- c(ptyOrd[i], ptyOrd[i] + r.res[i]*2, ptyOrd[i] + r.res[i]*2, ptyOrd[i]) - r.res[i]
      #       } else if (shapes[i] == "Diamond") {
      #         d <- 0.85 * r.res[i]*2*sqrt(2) # Diagonal of the square, since it's rotated. 0.8 lowers it because the vertices make it look bigger
      #         xi <- c(ptxOrd[i], ptxOrd[i] + d/2, ptxOrd[i] + d, ptxOrd[i] + d/2) - d/2
      #         yi <- c(ptyOrd[i], ptyOrd[i] + d/2, ptyOrd[i], ptyOrd[i] - d/2)
      #       } else if (shapes[i] == "Hexagon") {
      #         d <- 2*r.res[i]*2/sqrt(3) # Hexagon diagonal
      #         l <- 2*r.res[i]/sqrt(3) # Distance between two vertices
      #         p1 <- (d - l)/2 # y-axis movement step
      #         xi <- c(ptxOrd[i], ptxOrd[i] - r.res[i], ptxOrd[i] - r.res[i], ptxOrd[i], ptxOrd[i] + r.res[i], ptxOrd[i] + r.res[i])
      #         yi <- c(ptyOrd[i], ptyOrd[i] + p1, ptyOrd[i] + p1 + l, ptyOrd[i] + d, ptyOrd[i] + p1 + l, ptyOrd[i] + p1) - d/2
      #       } else { # Circles
      #         xi <- r.res[i] * sin(thetas) + ptxOrd[i]
      #         yi <- r.res[i] * cos(thetas) + ptyOrd[i]
      #       }
      #       
      #       polygon(xi, yi, col = cores[fills][i], lwd = bwd,
      #               border = bordersCol[i])
    }
    
    # Fill styles
    
    circsFills <- c(input$fill1, input$fill2, input$fill3, input$fill4, input$fill5)
    
    styles <- circsFills[fills]
    
    styles.angs <- c(n = -1, h = 0, v = 90, `d/` = 45, `d\\` = -45)
    
    if (any(styles != "n")) {
      nFills <- as.numeric(c(input$nFills1, input$nFills2, input$nFills3, input$nFills4, input$nFills5))
      nFills <- nFills[fills]
      fillCol <- c(input$fillCol1, input$fillCol2, input$fillCol3, input$fillCol4, input$fillCol5)
      fillCol <- fillCol[fills]
      angs.pattern <- styles.angs[match(styles, names(styles.angs))]
      for (i in nMin) {
        if (!is.na(styles[i]) && styles[i] != "n" && shapes[i] == "Circle") {
          #           ang <- if(styles[i] == "h") {
          #             0
          #           } else if (styles[i] == "v") {
          #             90
          #           } else if (styles[i] == "d/") {
          #             45
          #           } else if (styles[i] == "d\\") {
          #             -45
          #           }
          circStripes(r = r.res[i], n = nFills[i],
                      x = ptxOrd[i], y = ptyOrd[i],
                      angle = angs.pattern[i], col = fillCol[i])
        }
      }
      
      # Adding the cicle borders again...
      
      #       for (i in nMin) {
      #         xi <- r.res[i] * sin(thetas) + ptxOrd[i]
      #         yi <- r.res[i] * cos(thetas) + ptyOrd[i]
      #         polygon(xi, yi, col = NA, lwd = bordersWd[i],
      #                 border = bordersCol[i])
      #       }
    }
    
    # Using symbols is dangerous. The drawings are relative to the axis and drawings may be hard to position
    #     symbols(ptxOrd, ptyOrd, circles = r.res, inches = FALSE, 
    #             add = TRUE, bg = cores[fills], 
    #             fg = circBorder, lwd = input$circBorderWd) #Plotando os circulos
    
    # Making labels
    labType <- as.numeric(input$labType)
    
    resLab <-  if (labType == 0) {
      ""
    } else if (labType == 1) {
      res
    } else if (labType == 2) {
      amin3[match(res, amin)]
    } else if (labType == 3) {
      paste0(res, seq_len(nres) + input$labOffPos)
    }
    
    labOffX <- input$labOffX
    labOffY <- input$labOffY
    
    numLabs <- c(input$labCol1, input$labCol2, input$labCol3, input$labCol4, input$labCol5)
    
    text(ptxOrd + r.res*(sign(labOffX)) + labOffX,
         ptyOrd + r.res*(sign(labOffY)) + labOffY,
         label = resLab,
         cex = input$labCex/10, font = as.numeric(input$labFont),
         col = numLabs[fills])
    
    # Adding the residues number
    
    if (input$numShow == "Yes") {
      numOffY <- input$numOffY
      numOffX <- -input$numOffX
      
      numCols <- c(input$numCol1, input$numCol2, input$numCol3, input$numCol4, input$numCol5)
      
      text(ptxOrd - r.res*(sign(numOffX)) - numOffX, 
           ptyOrd - r.res*(sign(numOffY)) - numOffY, 
           labels = seq_len(nres) + input$numOffPos, cex = input$numCex/10, 
           font = as.numeric(input$numFont),
           col = numCols[fills]) 
    }
    
    if (input$showLeg == "Yes") {
      legLab <- c(input$leg1, input$leg2, input$leg3, input$leg4, input$leg5)
      legFill <- c(input$col1, input$col2, input$col3, input$col4, input$col5)
      if (!grepl("X", pepSeq)) {
        legLab <- legLab[1:4]
        legFill <- legFill[1:4]
      }
      
      # Defining a legend function to avoid repeating the same arguments below if necessary
      
      my.leg <- function(...) legend(x = r*input$legX, y = r*input$legY, 
                                     legend = legLab,
                                     bty = "n", yjust = 0.5, cex = as.numeric(input$legCex)/10, ...)
      
      # Now just calls with standard arguments and colouring fill
      
      my.leg(fill = legFill)
      
      # Adding pattern on legend
      
      angs.pattern <- styles.angs[match(circsFills, names(styles.angs))]
      #print(angs.pattern)
      if (any(!is.na(angs.pattern))) {
        
        legDen <- c(input$legDen1, input$legDen2, input$legDen3, input$legDen4, input$legDen5)
        
        my.leg(fill = ifelse(angs.pattern  == -1, adjustcolor("black", 0), "black"),
               density = ifelse(angs.pattern == -1, NA, legDen),
               angle = angs.pattern)
      }
    }
    
    # Adding title/name. Later to be on top
    
    if (input$showTitle == "Yes") {
      main <- input$txTitle
      cex.main <- input$cexTitle/10
      text(x = input$xTitle, y = input$yTitle, label = main,
           cex = cex.main, adj = c(0.5, 0.5), font = 2)
    }
    
    # Adding outer box to help exporting
    
    if (input$showBox == "Yes") {
      box("outer") 
    }
    
  }
  
  imgWd <- reactive(input$figWd)
  #getSeq <- reactive(input$seq)
  
  output$helicalPlot <- renderPlot({
    # The auto option is disbled by default to avoid unecessary plot creation 
    # when another projection is being tested, since this is rather slow
    # and require considerable resources when online
    if(input$autoWheel || is.null(lastWheel)) {
      helicalPlot()
      lastWheel <<- recordPlot()
    } else {
      #print(class(lastWheel))
      replayPlot(lastWheel)
    }
  },
  height = 600, width = imgWd)
  
  #shinyjs::disable("wheelsize")
  
  # Use sequence as title
  observeEvent(input$titleSeq, {
    updateSelectInput(session, "showTitle", selected = "Yes")
    updateTextInput(session, "txTitle", value = input$seq)
    nc <- nchar(input$seq)
    cexFit <- round(0.0325*nc^2 - 2.275*nc + 49.5, 1)
    updateNumericInput(session, "cexTitle", value = cexFit)
  })
  
  # Reseting and style buttons
  
  observeEvent(input$reset_input, {
    shinyjs::reset("side-panel")
  })
  
  observeEvent(input$reset_input2, {
    shinyjs::reset("settings")
  })
  
  
  observeEvent(input$style_gray, {
    shinyjs::updateColourInput(session = session, inputId = "col1", value = "black")
    updateTextInput(session = session, inputId = "circBorder1", value = "No")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol1", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol1", value = "#FFFFFF")
    shinyjs::updateColourInput(session = session, inputId = "col2", value = "gray30")
    updateTextInput(session = session, inputId = "circBorder1", value = "No")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol2", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol2", value = "#FFFFFF")
    shinyjs::updateColourInput(session = session, inputId = "col3", value = "gray50")
    updateTextInput(session = session, inputId = "circBorder3", value = "No")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol3", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol3", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "col4", value = "gray80")
    updateTextInput(session = session, inputId = "circBorder4", value = "No")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol4", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol4", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "col5", value = "white")
    updateTextInput(session = session, inputId = "circBorder5", value = "Yes")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol5", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol5", value = "#000000")
  })
  
  observeEvent(input$style_bw, {
    shinyjs::updateColourInput(session = session, inputId = "col1", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder1", value = "Yes")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol1", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol1", value = "#000000")
    updateSelectInput(session = session, inputId = "fill1", selected = "n")
    updateNumericInput(session = session, inputId = "nFills1", value = "4")
    shinyjs::updateColourInput(session = session, inputId = "fillCol1", value = "#FFFFFF")
    
    shinyjs::updateColourInput(session = session, inputId = "col2", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder1", value = "Yes")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol2", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol2", value = "#000000")
    updateSelectInput(session = session, inputId = "fill2", selected = "v")
    updateNumericInput(session = session, inputId = "nFills2", value = "4")
    shinyjs::updateColourInput(session = session, inputId = "fillCol2", value = "#000000")
    
    shinyjs::updateColourInput(session = session, inputId = "col3", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder3", value = "Yes")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol3", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol3", value = "#000000")
    updateSelectInput(session = session, inputId = "fill3", selected = "h")
    updateNumericInput(session = session, inputId = "nFills3", value = "4")
    shinyjs::updateColourInput(session = session, inputId = "fillCol3", value = "#000000")
    
    shinyjs::updateColourInput(session = session, inputId = "col4", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder4", value = "Yes")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol4", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol4", value = "#000000")
    updateSelectInput(session = session, inputId = "fill4", selected = "d/")
    updateNumericInput(session = session, inputId = "nFills4", value = "3")
    shinyjs::updateColourInput(session = session, inputId = "fillCol4", value = "#000000")
    
    shinyjs::updateColourInput(session = session, inputId = "col5", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder5", value = "Yes")
    shinyjs::updateColourInput(session = session, inputId = "circBorderCol5", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol5", value = "#000000")
    updateSelectInput(session = session, inputId = "fill5", selected = "d\\")
    updateNumericInput(session = session, inputId = "nFills5", value = "3")
    shinyjs::updateColourInput(session = session, inputId = "fillCol5", value = "#000000")
    
    updateNumericInput(session = session, inputId = "maxlinegray", value = 1)
    updateNumericInput(session = session, inputId = "minlinegray", value = 1)
    updateNumericInput(session = session, inputId = "conLineWd", value = 1)
  })
  
  observeEvent(input$style_labdown, {
    current <- input$labType
    
    if (current == "0") {
      updateSelectInput(session = session, inputId = "numShow", selected = "Yes")
      updateSelectInput(session = session, inputId = "labType", selected = "1")
      updateNumericInput(session = session, inputId = "labOffY", value = 0)
      updateNumericInput(session = session, inputId = "lnumOffY", value = -0.3)
    } else if (current == "1") {
      updateSelectInput(session = session, inputId = "labType", selected = "2")
    } else if (current == "2") {
      updateSelectInput(session = session, inputId = "numShow", selected = "No")
      updateSelectInput(session = session, inputId = "labType", selected = "3")
      updateNumericInput(session = session, inputId = "labOffY", value = -0.3)
    } else if (current == "3") {
      if (input$labOffY == -0.3) {
        updateNumericInput(session = session, inputId = "labOffY", value = 0)
      } else {
        updateSelectInput(session = session, inputId = "numShow", selected = "No")
      updateSelectInput(session = session, inputId = "labType", selected = "0")
      }
    } 
    shinyjs::updateColourInput(session = session, inputId = "labCol1", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol2", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol3", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol4", value = "#000000")
    shinyjs::updateColourInput(session = session, inputId = "labCol5", value = "#000000")
    
  })
  
  observeEvent(input$style_round, {
    if (input$innRatio != 0) {
      updateNumericInput(session, "innRatio", value = 0)
    } else {
      updateNumericInput(session, "innRatio", value = 0.4)  
    }
  })
  
  observeEvent(input$style_legend, {
    if (input$showLeg == "No") {
      updateSelectInput(session, "showLeg", selected = "Yes")
      updateNumericInput(session, "figWd", value = 775)
      updateNumericInput(session, "marR", value = 12)
      updateNumericInput(session, "marT", value = 0)
    } else {
      updateSelectInput(session, "showLeg", selected = "No")
      updateNumericInput(session, "figWd", value = 600)
      updateNumericInput(session, "marR", value = 0)
      updateNumericInput(session, "marT", value = 0)
    }
  })
  
  observeEvent(input$style_circle, {
    updateSelectInput(session = session, inputId = "shp1", selected = "Circle")
    updateSelectInput(session = session, inputId = "shp2", selected = "Circle")
    updateSelectInput(session = session, inputId = "shp3", selected = "Circle")
    updateSelectInput(session = session, inputId = "shp4", selected = "Circle")
    updateSelectInput(session = session, inputId = "shp5", selected = "Circle")
  })
  
  observeEvent(input$style_squares, {
    updateSelectInput(session = session, inputId = "shp1", selected = "Square")
    updateSelectInput(session = session, inputId = "shp2", selected = "Square")
    updateSelectInput(session = session, inputId = "shp3", selected = "Square")
    updateSelectInput(session = session, inputId = "shp4", selected = "Square")
    updateSelectInput(session = session, inputId = "shp5", selected = "Square")
  })
  
  observeEvent(input$style_shapes, {
    updateSelectInput(session = session, inputId = "shp1", selected = "Circle")
    updateSelectInput(session = session, inputId = "shp2", selected = "Square")
    updateSelectInput(session = session, inputId = "shp3", selected = "Triangle")
    updateSelectInput(session = session, inputId = "shp4", selected = "Diamond")
    updateSelectInput(session = session, inputId = "shp5", selected = "Hexagon")
  })
  
  observeEvent(input$bt_rdata, {
    save.image("Image.RData")
  })
  
  observeEvent(input$bt_rdataload, {
    load("Image.RData")
  })
  
  observeEvent(input$bt_force, {
    updateSelectInput(session, "circBorder1", selected = "Yes")
  })
 
  
  output$bt_export <- downloadHandler(
    filename = function() {
      "HelicalWheelStyle.hws"
    },
    content = function(file) {
      vlist <- reactiveValuesToList(input)
      # Remove this values because they aren't inputs and can be dangerous
      vlist[c("fileImport", "shinyjs-resettable-settings")] <- NULL
      #print(vlist)
      inputsList <- names(vlist)
      exportVars <- paste0(inputsList, "=", sapply(inputsList, function(inpt) input[[inpt]]))
      write(exportVars, file)
    })
  
  importFile <- reactive({
    
    inFile <- input$fileImport
    #print("a")
    if (is.null(inFile))
      return(NULL)
    
    lines <- readLines(inFile$datapath)
    out <- lapply(lines, function(l) unlist(strsplit(l, "=")))
    #shinyjs::reset("fileImport")
    return(out)
  })
  
  observe({
    imp <- importFile()

    if (!is.null(imp)) {
      imp <- imp[order(sapply(imp, function(o) nchar(o[1])))]

      for (inpt in imp) session$sendInputMessage(inpt[1], list(value = inpt[2]))
      
      # This is a workaround to make sure reuploading the file will work properly. 
      Sys.sleep(1)
      session$sendCustomMessage(type = "resetFileInputHandler", 'fileImport')
    }
  })
  
  fname <- "HelicalWheelProjection"
  
  output$plotDown <- downloadHandler(
    filename = function() {
      paste(fname, tolower(input$expFormat), sep = ".")
    },
    content = function(file) {
      dpi <- input$expDPI
      dims <- as.numeric(dpi)*c(input$figWd, 600)/72
      if (input$expFormat == "PNG") {
        png(file, res = dpi, width = dims[1], height = dims[2])
      } else if (input$expFormat == "TIFF") {
        tiff(file, compression = "lzw", res = dpi, width = dims[1], height = dims[2])
      } else if (input$expFormat == "PDF") {
        pdf(file, width = input$figWd/72, height = 600/72)
      } else if (input$expFormat == "JPG") {
        jpeg(file, res = dpi, width = dims[1], height = dims[2], quality = 100)
      }
      par(mar = c(0, 0, input$marT, input$marR)+0.1, xpd = TRUE)
      helicalPlot()
      dev.off()
    }
  )
  
  fnameNet <- "HelicalNetProjection"
  
  output$plotDownNet <- downloadHandler(
    filename = function() {
      paste(fnameNet, tolower(input$expFormatNet), sep = ".")
    },
    content = function(file) {
      dpi <- input$expDPINet
      dims <- as.numeric(dpi)*c(300, getHeight())/72
      ImgFormat <- input$expFormatNet 
      if (ImgFormat == "PNG") {
        png(file, res = dpi, width = dims[1], height = dims[2])
      } else if (ImgFormat == "TIFF") {
        tiff(file, compression = "lzw", res = dpi, width = dims[1], height = dims[2])
      } else if (ImgFormat == "PDF") {
        pdf(file, width = 300/72, height = getHeight()/72)
      } else if (ImgFormat == "JPG") {
        jpeg(file, res = dpi, width = dims[1], height = dims[2], quality = 100)
      }
      par(mar=c(0,0,0,0) + 0.1)
      netPlot()
      dev.off()
    }
  )
})
