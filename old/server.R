library(shiny)

# These objects aren't reactive, they are written here so they are only created once for the whole app, and not for every session.

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

#   rotvec <- function (v, n) {
#     # Taken from wavethresh::guyrot
#     l <- length(v)
#     n <- n%%l
#     if (n == 0) 
#       return(v)
#     tmp <- v[(l - n + 1):l]
#     v[(n + 1):l] <- v[1:(l - n)]
#     v[1:n] <- tmp
#     v
#   }

shinyServer(function(input, output, session) {
  
  ## All aminoacids and their classification
  #amin <- c("R", "H", "K", "D", "E", "S", "T", "N", "Q", "C", "U", "P", "A", "G", "V", "I", "L", "M", "F", "Y", "W", "X")
  #hb   <- c( 1 ,  1 ,  1 ,  2 ,  2 ,  3 ,  3 ,  3 ,  3 ,  3 ,  3 ,  4 ,  4 ,  4 ,  4 ,  4 ,  4 ,  4 ,  4 ,  4 ,  4  , 5 )
  ##        charged/basic    -acid-    ------------polar/uncharged---------    -------------hydrophibic-------------
  # The values above are the default ones.
  # Below is the code to load them from the inputs and give them their proper codes
  # The reactive expressions are evaluated inside netPlot() and wheelPlot() to create the amin and hb variables
  getAmin <- reactive(toupper(unlist(strsplit(paste0(c(input$grp1, input$grp2, input$grp3, input$grp4, input$grp5), collapse = ""), ""))))
  
  gethb <- reactive(c(rep(1, nchar(input$grp1)), rep(2, nchar(input$grp2)), rep(3, nchar(input$grp3)), 
          rep(4, nchar(input$grp4)), rep(5, nchar(input$grp5))))
  
  
  # The 'X' aminoacid represents any aminoacid. All characters not present in amin will be replaced with X, unless a 3-letter format is detected
  
  amin3 <-c(R = "Arg", H = "His", K = "Lys", D = "Asp", E = "Glu", S = "Ser", T = "Thr", N = "Asn",
            Q = "Gln", C = "Cys", U = "Sec", G = "Gly", P = "Pro", A = "Ala", V = "Val", I = "Ile",
            L = "Leu", M = "Met", F = "Phe", Y = "Tyr", W = "Trp")
  
  # Peptide groups in terms of bonds, some do miss
  
#   nonpolar <- c("V", "I", "L", "M", "F", "Y", "W")
#   acid <- c("D", "E")
#   basic <- c("R", "H", "K")
#   hydrobond <- c("S", "N", "Q")
  
  getBonds <- reactive(lapply(list(input$grpNonpolar,
                            input$grpAcid,
                            input$grpBasic,
                            input$grpHydro), function(i) unlist(strsplit(i, ""))))
  
  #### Code for Helical Net ####
  
  getHeight <- function() {
    pepSeq <- input$seq
    n <- nchar(pepSeq)
    rws <- (n*1.5)/5.4
    # The 120 factor is empirical, it creates the final height 
    # to make sure the symbols have proportional dimensions
    height <- 120*rws
    height
  }
  
  getWidth <- function() {
    300 + input$netWidth  
  }
  
  lastNet <- NULL
  
  countFunc <- reactive({
    seq <- input$seq
    output$resCount <- renderUI(helpText(seq, style = "font-family:monospace; margin: 0 12px; font-size: 16px; "))
    n <- nchar(seq)
    counts <- paste0(paste(rep(c(1:9, 0), 10)[(1:n) + (input$numOffPos %% 10)], collapse = ""), "(", n, ")")
    
    output$seqMono <- renderUI(helpText(counts,
                                        style = "font-family:monospace; margin: 0 12px; font-size: 16px; "))
  })

  netPlot <- function() {
    pepSeq <- toupper(input$seq)
    amin <- getAmin()
    hb <- gethb()
    n <- nchar(pepSeq)
    
    bondsInpt <- getBonds()
    
    nonpolar <- bondsInpt[[1]]
    acid <- bondsInpt[[2]]
    basic <- bondsInpt[[3]]
    hydrobond <- bondsInpt[[4]]
    
    #output$resNumber <- renderUI(helpText(paste("Number of residues:", n)))
    #resCount <- paste(rep(c(1:9, 0), 10)[1:n], collapse = "")
    
    countFunc()
    
    intercept <- 0
    
    ###### Projection properties inputs
    prop_fact <- input$netProp #3
    #diameter <- input$netDiameter #4.6
    #pitch <- input$netPitch/prop_fact # 5.4/prop_fact
    #trans <- input$netTrans/prop_fact # 1.5/prop_fact
    y0 = input$netStartOff/prop_fact #0
    perturn <- input$netPerTurn #3.6
    netDir <- as.numeric(input$netDirection)
    #######
    
    ####### Other projection graphic inputs
    marl <- input$netPadL 
    marr <- input$netPadR
    padtop <- input$netPadTop
    padbot <- input$netPadBot
    netwd <- input$netWidth/60
    #######
    
    ####### This values are in Angstrons for the alpa helix, but in practice they can be fixed
    # The residues per turn is enough to actually change the position of the residues
    diameter <- 4.6
    pitch <- 5.4/prop_fact
    trans <- 1.5/prop_fact
    #######
    
    x0 <- diameter/5
    
    ymax <- max(n)*trans+y0
    ny <- seq(from = ymax, by = -trans, length.out = n)
    nx1 <- (ny[1]-intercept)*(diameter/pitch)
    nx <- seq(from = nx1, by = -netDir*diameter/perturn, length.out = n)
    #nx <- (ny-intercept)*(diameter/pitch)
    nx <- ifelse(nx > diameter, nx-(nx%/%diameter)*diameter, nx)
    
    par(mar=c(0,0,0,0) + 0.1)
#    plot(1, type = "n", xlim = c(0, 4.6), ylim = c(min(ny)-trans, ymax+trans))
    plot(1, type = "n", xlim = c(-marl, diameter + marr + netwd), 
         ylim = c(min(ny)-trans-padbot, ymax+trans+padtop),
         xaxt = "n", yaxt = "n", xlab = "", ylab = "",
         xaxs = "i", yaxs = "i", frame.plot = input$showBoxNet == "Yes")

    xymat <- do.call(rbind, lapply(seq_len(n), function(i) {
      if (diameter - nx[i] < x0 && diameter - nx[i] > 0) {
        # Close to right border
        matrix(c(nx[i]-diameter, nx[i], rep(ny[i], 2), rep(i, 2), c(0, 1)),
               ncol = 4)
      } else if(nx[i]-intercept < x0) {
        # Close to left border
        matrix(c(nx[i], nx[i]+diameter, rep(ny[i], 2), rep(i, 2), c(1, 0)),
               ncol = 4)
      } else {
        matrix(c(nx[i], ny[i], i, 1), ncol = 4)
      }
    }))
    
    ptx <- xymat[,1]
    pty <- xymat[,2]
    ptn <- xymat[,3]
    ptr <- xymat[,4]
    
    # Setting the sizes of the polygons. 
    # The /2/2 divisions just allows for the same input range even with different output sizes for net and wheel
    lfrac <- input$circsize/2
    l <- diameter*lfrac/2
    # One prop value is created for each unique residue, but it's repeated for residues which will be plotted twice using [ptn]
    l.prop <- sqrt(seq(input$circprop, 1, length.out = n))[ptn]
    
    l.prop <- l.prop/max(l.prop)
    
    # Adding guides behind
    
    if (input$showNetGuide == "Yes") {
      ltyNetGuide <- as.numeric(input$netGuideLty)
      lwdNetGuide <- as.numeric(input$netGuideLwd)
      colNetGuide <- gray(1-as.numeric(input$netGuideCol))
      
      #ord <- order(pty, sort(ptx))
      lines(ptx[ptr==1], pty[ptr==1], lty = ltyNetGuide, lwd = lwdNetGuide, col = colNetGuide) 
    }
    
    res <- unlist(strsplit(pepSeq, ""))[ptn]
    fills <- hb[match(res, amin)]
    
    shp <- c(input$shp1, input$shp2, input$shp3, input$shp4, input$shp5)
    cores <- c(input$col1, input$col2, input$col3, input$col4, input$col5)
    bordersShow <- c(input$circBorder1, input$circBorder2,
                     input$circBorder3, input$circBorder4, input$circBorder5)
    
    if (any(bordersShow == "Yes")) {
      bordersCol <- c(input$circBorderCol1, input$circBorderCol2,
                      input$circBorderCol3, input$circBorderCol4, input$circBorderCol5)
      bordersCol[bordersShow != "Yes"] <- NA
      bordersCol <- bordersCol[fills]
      bordersWd <- as.numeric(c(input$circBorderWd1, input$circBorderWd2,
                                input$circBorderWd3, input$circBorderWd4, input$circBorderWd5))
      bordersCol[bordersWd <= 0] <- NA
      bordersWd <- ifelse(bordersWd == 0, 1, bordersWd)
      bordersWd <- bordersWd[fills]
    } else {
      bordersCol <- rep(NA, length(ptn))
      bordersWd <- rep(1, length(ptn))
    }
    
    shapes <- shp[fills]
    
    nres <- length(ptn)
    numLabs <- c(input$labCol1, input$labCol2, input$labCol3, input$labCol4, input$labCol5)
    
    # Making labels
    labType <- as.numeric(input$labType)
    
    resLab <-  if (labType == 0) {
      ""
    } else if (labType == 1) {
      res
    } else if (labType == 2) {
      amin3[res]
    } else if (labType == 3) {
      paste0(res, ptn+input$numOffPos)
    } else if (labType == 4) {
      ptn+input$numOffPos
    }
    
    bondWd <- c(nonpolar = input$bond1Wd, acba = input$bond2Wd, hydro = input$bond3Wd)
    bondTy <- c(nonpolar = input$bond1Ty, acba = input$bond2Ty, hydro = input$bond3Ty)
    bondCol <- c(nonpolar = input$bond1Col, acba = input$bond2Col, hydro = input$bond3Col)
    
    #Looping over residues to add their polygons/bonds
    
    for (i in seq_len(nres)) {
      # Detecting Bonds 
      
      if (input$netShowInteractions == "Yes") {
      
        bond3 <- TRUE
        next3 <- which(ptn == ptn[i] + 3)[1]
        
        if (res[i] %in% nonpolar && res[next3] %in% nonpolar) {
          bond3type <- "nonpolar"
        } else if ((res[i] %in% basic && res[next3] %in% acid) ||
                   (res[i] %in% acid && res[next3] %in% basic)) {
          bond3type <- "acba"
        } else if (res[i] %in% hydrobond && res[next3] %in% hydrobond) {
          bond3type <- "hydro"
        } else {
          bond3 <- FALSE
          bond3type <- "none"
        } 
        
        bond4 <- TRUE
        next4 <- which(ptn == ptn[i] + 4)[1]
        if (res[i] %in% nonpolar && res[next4] %in% nonpolar) {
          bond4type <- "nonpolar"
        } else if ((res[i] %in% basic && res[next4] %in% acid) ||
                   (res[i] %in% acid && res[next4] %in% basic)) {
          bond4type <- "acba"
        } else if (res[i] %in% hydrobond && res[next4] %in% hydrobond) {
          bond4type <- "hydro"
        } else {
          bond4 <- FALSE
          bond4type <- "none"
        } 
        # Adding bonds
       
        if (input$netDirection == "1") {
          # The min/max choice for broken bonds depends on the direction of the net
          # We select the min or max functions to be the extremes depending on the direction
          ext1 <- max
          ext2 <- min
        } else {
          ext1 <- min
          ext2 <- max
        }
        
        lty3 <- if(bond3 && nchar(bondTy[bond3type]) > 1) bondTy[bond3type] else as.numeric(bondTy[bond3type])
        lty4 <- if(bond4 && nchar(bondTy[bond4type]) > 1) bondTy[bond4type] else as.numeric(bondTy[bond4type])
        
        if (bond3) {
          segments(x0 = ptx[i], y0 = pty[i]-l,
                   x1 = ptx[i]-((diameter/perturn)*3-diameter)*netDir,
                   y1 = pty[i]-3*trans+l,
                   lwd = as.numeric(bondWd[bond3type]), lty = lty3,
                   col = bondCol[bond3type])
        }
  
        if (bond4) {
          segments(x0 = ptx[i], y0 = pty[i]-l,
                   x1 = ptx[i]-((diameter/perturn)*4-diameter)*netDir,
                   y1 = pty[i]-4*trans+l,
                   lwd = as.numeric(bondWd[bond4type]), lty = lty4,
                   col = bondCol[bond4type])
        }
      }
        
      # Adding Residues Polygons
      bwd <- bordersWd[i]
      addPolygon(shape = shapes[i], 
                 x = ptx[i], y = pty[i],
                 size = l*l.prop[i], col = cores[fills][i],
                 lwd = bwd, border = bordersCol[i])
      
    }
    
    # Fill styles
    
    circsFills <- c(input$fill1, input$fill2, input$fill3, input$fill4, input$fill5)
    
    styles <- circsFills[fills]
    
    styles.angs <- c(n = -1, h = 0, v = 90, `d/` = 45, `d\\` = -45)
    
    if (all(!is.na(styles)) && any(styles != "n")) {
      nFills <- as.numeric(c(input$nFills1, input$nFills2, input$nFills3, input$nFills4, input$nFills5))
      nFills <- nFills[fills]
      fillCol <- c(input$fillCol1, input$fillCol2, input$fillCol3, input$fillCol4, input$fillCol5)
      fillCol <- fillCol[fills]
      angs.pattern <- styles.angs[match(styles, names(styles.angs))]
      for (i in seq_len(nres)) {
        if (!is.na(styles[i]) && styles[i] != "n" && shapes[i] == "Circle") {
          circStripes(r = l*l.prop[i], n = nFills[i],
                      x = ptx[i], y = pty[i],
                      angle = angs.pattern[i], col = fillCol[i])
        }
      }
    }
    
    labOffX <- input$labOffX/5
    labOffY <- input$labOffY/5
    
    text(ptx + l*(sign(labOffX)) + labOffX, 
         pty + l*(sign(labOffY)) + labOffY,
         label = resLab,
         cex = input$labCex/10, font = as.numeric(input$labFont),
         col = numLabs[fills])
    
    # Adding number labels
    
    if (input$numShow == "Net" || input$numShow == "Both") {
      numOffY <- input$numOffY
      numOffX <- -input$numOffX
      
      numCols <- c(input$numCol1, input$numCol2, input$numCol3, input$numCol4, input$numCol5)
      
      text(ptx - l*(sign(numOffX))*l.prop - numOffX, 
           pty - l*(sign(numOffY))*l.prop - numOffY, 
           labels = ptn + input$numOffPos, cex = input$numCex/10, 
           font = as.numeric(input$numFont),
           col = numCols[fills]) 
    }
    
    ### Adding a white rectangular to cover residues plotted outside of the region of interest
    
    polygon(x = c(diameter + marr, diameter + marr + netwd - 0.05, diameter + marr + netwd - 0.05, diameter + marr),
            y = c(ymax + trans + padtop - 0.05, ymax + trans + padtop - 0.05, min(ny) - trans - padbot, min(ny) - trans-padbot),
            col = "white", border = NA)
    
    ### Adding vertical cylinder lines
    
    if (input$netShowLimits == "Yes") abline(v=c(intercept, diameter), lty = 2)
    
    ### Legend adding code:
    
    # The 2 lines below check which characters given as a Sequence are present in each groups, in order to only add the necessary legends
    grpSplit <- split(getAmin(), gethb())
    groupsPresent <- sapply(grpSplit, function(gp) any(strsplit(pepSeq, "")[[1]] %in% gp))
    
    if (input$showLegNet == "Yes") {
      legLab <- c(input$leg1, input$leg2, input$leg3, input$leg4, input$leg5)[groupsPresent]
      legFill <- c(input$col1, input$col2, input$col3, input$col4, input$col5)[groupsPresent]
      #       if (!grepl("X", pepSeq)) {
      #         legLab <- legLab[1:4]
      #         legFill <- legFill[1:4]
      #       }
      
      # Defining a legend function to avoid repeating the same arguments below if necessary
      
      my.leg <- function(...) legend(x = diameter + input$legXNet, 
                                     y = ymax+trans+padtop + input$legYNet, 
                                     legend = legLab,
                                     bty = "n", yjust = 0.5, cex = as.numeric(input$legCex)/10, ...)
      
      # Now just calls with standard arguments and colouring fill
      
      #polygon(x = input$legX, y = input$legY)
      
      my.leg(fill = legFill)
      
      # Adding pattern on legend
      
      angs.pattern <- styles.angs[match(circsFills, names(styles.angs))]
      
      if (any(!is.na(angs.pattern))) {
        
        legDen <- c(input$legDen1, input$legDen2, input$legDen3, input$legDen4, input$legDen5)
        
        my.leg(fill = ifelse(angs.pattern  == -1, adjustcolor("black", 0), "black"),
               density = ifelse(angs.pattern == -1, NA, legDen),
               angle = angs.pattern)
      }
    }
    
    ### End of legend code.
    
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
      if(input$autoNet || is.null(lastNet)) {
        netPlot()
        lastNet <<- TRUE
      } else {
        plot(1:10, 1:10, type = "n", xaxt = "n", yaxt = "n", frame.plot = FALSE,
             xlab = "", ylab = "", main = NA)
        text(5, 7, cex=1.2,
             labels = "Check the 'Automatic Preview'\nbox to see the results.\n\nThis can be a bit slow.")
      }
    },
    height = getHeight, width = getWidth)
  
  
  #### Code for Helical Wheel ####
  
  # Wheel radius. Changing this just messes with some relative positions offsets.
  r <- 10 #input$wheelsize # Raio do circulo grande  da figura
  
  lastWheel <- NULL
  
  helicalPlot <- function() {
    amin <- getAmin()
    hb <- gethb()
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
    
    #output$resNumber <- renderUI(helpText(paste("Number of residues:", nres)))
    #resCount <- paste(rep(c(1:9, 0), 10)[1:nres], collapse = "")
    
    countFunc()
    
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
    
    r.res.prop <- sqrt(r.res.prop)
    # The circles should have the area, not the radius, proportional do the position. This fixes it based on A = pi * r * r
    # It's based on circles, but the sqrt idea is the same for other polygons.

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
    
    bordersShow <- c(input$circBorder1, input$circBorder2,
                     input$circBorder3, input$circBorder4, input$circBorder5)
    
    if (any(bordersShow == "Yes")) {
      bordersCol <- c(input$circBorderCol1, input$circBorderCol2,
                      input$circBorderCol3, input$circBorderCol4, input$circBorderCol5)
      bordersCol[bordersShow != "Yes"] <- NA
      bordersCol <- bordersCol[fills]
      bordersWd <- as.numeric(c(input$circBorderWd1, input$circBorderWd2,
                                input$circBorderWd3, input$circBorderWd4, input$circBorderWd5))
      bordersCol[bordersWd <= 0] <- NA
      bordersWd <- ifelse(bordersWd == 0, 1, bordersWd)
      bordersWd <- bordersWd[fills]
    } else {
      bordersCol <- rep(NA, length(nMin))
      bordersWd <- rep(1, length(nMin))
    }
    
    # Ploting the polygons

    shp <- c(input$shp1, input$shp2, input$shp3, input$shp4, input$shp5)
    
    shapes <- shp[fills]
    
    for (i in 1:nres) {
      bwd <- bordersWd[i]
      #if (is.na(bwd)) next # Avoiding errors with invalid lwd. NA throwns an error, skip because there's no circle to plot
      bordersWd <- ifelse(bordersWd == 0, 1, bordersWd)
      
      addPolygon(shape = shapes[i], 
                 x = ptxOrd[i], y = ptyOrd[i],
                 size = r.res[i], col = cores[fills][i],
                 lwd = bwd, border = bordersCol[i])
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

          circStripes(r = r.res[i], n = nFills[i],
                      x = ptxOrd[i], y = ptyOrd[i],
                      angle = angs.pattern[i], col = fillCol[i])
        }
      }
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
      amin3[res]
    } else if (labType == 3) {
      paste0(res, seq_len(nres) + input$numOffPos)
    } else if (labType == 4) {
      seq_len(nres)
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
    
    if (input$numShow == "Wheel" || input$numShow == "Both") {
      numOffY <- input$numOffY
      numOffX <- -input$numOffX
      
      numCols <- c(input$numCol1, input$numCol2, input$numCol3, input$numCol4, input$numCol5)
      
      text(ptxOrd - r.res*(sign(numOffX)) - numOffX, 
           ptyOrd - r.res*(sign(numOffY)) - numOffY, 
           labels = seq_len(nres) + input$numOffPos, cex = input$numCex/10, 
           font = as.numeric(input$numFont),
           col = numCols[fills]) 
    }
    
    
    ### Legend adding code:
    
    # The 2 lines below check which characters given as a Sequence are present in each groups, in order to only add the necessary legends
    grpSplit <- split(getAmin(), gethb())
    groupsPresent <- sapply(grpSplit, function(gp) any(strsplit(pepSeq, "")[[1]] %in% gp))

    if (input$showLeg == "Yes") {
      legLab <- c(input$leg1, input$leg2, input$leg3, input$leg4, input$leg5)[groupsPresent]
      legFill <- c(input$col1, input$col2, input$col3, input$col4, input$col5)[groupsPresent]
#       if (!grepl("X", pepSeq)) {
#         legLab <- legLab[1:4]
#         legFill <- legFill[1:4]
#       }
      
      # Defining a legend function to avoid repeating the same arguments below if necessary
      
      my.leg <- function(...) legend(x = r*input$legX, y = r*input$legY, 
                                     legend = legLab,
                                     bty = "n", yjust = 0.5, cex = as.numeric(input$legCex)/10, ...)
      
      # Now just calls with standard arguments and colouring fill

      my.leg(fill = legFill)
      
      # Adding pattern on legend
      
      angs.pattern <- styles.angs[match(circsFills, names(styles.angs))]
      
      if (any(!is.na(angs.pattern))) {
        
        legDen <- c(input$legDen1, input$legDen2, input$legDen3, input$legDen4, input$legDen5)
        
        my.leg(fill = ifelse(angs.pattern  == -1, adjustcolor("black", 0), "black"),
               density = ifelse(angs.pattern == -1, NA, legDen),
               angle = angs.pattern)
      }
    }
    
    ### End of legend code.
    
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
    # The auto option is disabled by default to avoid unnecessary plot creation 
    # when another projection is being tested, since this is rather slow
    # and require considerable resources when online
    if(input$autoWheel || is.null(lastWheel)) {
      helicalPlot()
      lastWheel <<- recordPlot()
    } else {
      plot(1:10, 1:10, type = "n", xaxt = "n", yaxt = "n", frame.plot = FALSE,
           xlab = "", ylab = "", main = NA)
      text(5, 7, cex = 1.5,
           labels = "Check the 'Automatic Preview'\nbox to see the results.\n\nThis can be a bit slow.")
    }
  },
  height = 600, width = imgWd)
  
  #shinyjs::disable("period")
  
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
  
  observeEvent(input$grpReset, {
    updateTextInput(session, "grp1", value = "RHK")
    updateTextInput(session, "grp2", value = "DE")
    updateTextInput(session, "grp3", value = "STNQC")
    updateTextInput(session, "grp4", value = "AGVILMFYWP")
    updateTextInput(session, "grp5", value = "X")
  })
  
  observeEvent(input$grpBondReset, {
    updateTextInput(session, "grpNonpolar", value = "VILMFYW")
    updateTextInput(session, "grpHydro", value = "STNQY")
    updateTextInput(session, "grpAcid", value = "DE")
    updateTextInput(session, "grpBasic", value = "RHK")
  })
  
  
  
  
  grpUpdate <- function(newGrps) {
      for (i in 1:5) {
        colourpicker::updateColourInput(session, paste0("numCol", i), label = newGrps[i])
        updateTextInput(session, paste0("leg", i), label = newGrps[i])
        updateNumericInput(session, paste0("legDen", i), label = newGrps[i])
        updateSelectInput(session, paste0("grp", i), label = newGrps[i])
        updateSelectInput(session, paste0("shp", i), label = newGrps[i])
        colourpicker::updateColourInput(session, paste0("col", i), label = newGrps[i])
        updateSelectInput(session, paste0("circBorder", i), label = newGrps[i])
        colourpicker::updateColourInput(session, paste0("labCol", i), label = newGrps[i])
        updateSelectInput(session, paste0("fill", i), label = newGrps[i])
      }
  }
  
  observeEvent(input$grpUpdUi, {
    grpUpdate(c(input$grp1Lab, input$grp2Lab, input$grp3Lab, input$grp4Lab, input$grp5Lab))
  })
  
  observeEvent(input$grpLabDefault, {
    defaults <- c("Polar / Basic", "Polar / Acidic", "Polar / Uncharged", "Nonpolar", "Unkown Residue")
    for (i in 1:5) updateTextInput(session, paste0("grp", i, "Lab"), value = defaults[i])
    grpUpdate(defaults)
  })
  
  observeEvent(input$grpLabToLeg, {
    for (i in 1:5) updateTextInput(session, paste0("leg", i), value = input[[paste0("grp", i, "Lab")]])
  })
  
  
  observeEvent(input$netAutoMar, {
    nc <- nchar(input$seq)
    updateNumericInput(session, "yTitleNet", value = round(nc*0.51+0.5, 1))
    updateNumericInput(session, "netPadTop", value = 0.5)
  })
  
  observeEvent(input$style_gray, {
    colourpicker::updateColourInput(session = session, inputId = "col1", value = "black")
    updateTextInput(session = session, inputId = "circBorder1", value = "No")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol1", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol1", value = "#FFFFFF")
    colourpicker::updateColourInput(session = session, inputId = "col2", value = "gray30")
    updateTextInput(session = session, inputId = "circBorder2", value = "No")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol2", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol2", value = "#FFFFFF")
    colourpicker::updateColourInput(session = session, inputId = "col3", value = "gray50")
    updateTextInput(session = session, inputId = "circBorder3", value = "No")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol3", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol3", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "col4", value = "gray80")
    updateTextInput(session = session, inputId = "circBorder4", value = "No")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol4", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol4", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "col5", value = "white")
    updateTextInput(session = session, inputId = "circBorder5", value = "Yes")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol5", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol5", value = "#000000")
  })
  
  terColors <- substr(terrain.colors(5), 0, 7)
  observeEvent(input$style_terrain, {
    colourpicker::updateColourInput(session = session, inputId = "col1", value = terColors[1])
    updateTextInput(session = session, inputId = "circBorder1", value = "No")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol1", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol1", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "col2", value = terColors[2])
    updateTextInput(session = session, inputId = "circBorder2", value = "No")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol2", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol2", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "col3", value = terColors[3])
    updateTextInput(session = session, inputId = "circBorder3", value = "No")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol3", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol3", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "col4", value = terColors[4])
    updateTextInput(session = session, inputId = "circBorder4", value = "No")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol4", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol4", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "col5", value = terColors[5])
    updateTextInput(session = session, inputId = "circBorder5", value = "Yes")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol5", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol5", value = "#000000")
  })
  
  observeEvent(input$style_bw, {
    colourpicker::updateColourInput(session = session, inputId = "col1", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder1", value = "Yes")
    updateNumericInput(session = session, inputId = "circBorderWd1", value = 5)
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol1", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol1", value = "#000000")
    updateSelectInput(session = session, inputId = "fill1", selected = "n")
    updateNumericInput(session = session, inputId = "nFills1", value = "4")
    colourpicker::updateColourInput(session = session, inputId = "fillCol1", value = "#FFFFFF")
    
    colourpicker::updateColourInput(session = session, inputId = "col2", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder2", value = "Yes")
    updateNumericInput(session = session, inputId = "circBorderWd2", value = 4)
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol2", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol2", value = "#000000")
    updateSelectInput(session = session, inputId = "fill2", selected = "v")
    updateNumericInput(session = session, inputId = "nFills2", value = "4")
    colourpicker::updateColourInput(session = session, inputId = "fillCol2", value = "#000000")
    
    colourpicker::updateColourInput(session = session, inputId = "col3", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder3", value = "Yes")
    updateNumericInput(session = session, inputId = "circBorderWd3", value = 3)
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol3", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol3", value = "#000000")
    updateSelectInput(session = session, inputId = "fill3", selected = "h")
    updateNumericInput(session = session, inputId = "nFills3", value = "5")
    colourpicker::updateColourInput(session = session, inputId = "fillCol3", value = "#000000")
    
    colourpicker::updateColourInput(session = session, inputId = "col4", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder4", value = "Yes")
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol4", value = "#000000")
    updateNumericInput(session = session, inputId = "circBorderWd4", value = 2)
    colourpicker::updateColourInput(session = session, inputId = "labCol4", value = "#000000")
    updateSelectInput(session = session, inputId = "fill4", selected = "d/")
    updateNumericInput(session = session, inputId = "nFills4", value = "6")
    colourpicker::updateColourInput(session = session, inputId = "fillCol4", value = "#000000")
    
    colourpicker::updateColourInput(session = session, inputId = "col5", value = "#FFFFFF")
    updateTextInput(session = session, inputId = "circBorder5", value = "Yes")
    updateNumericInput(session = session, inputId = "circBorderWd5", value = 1)
    colourpicker::updateColourInput(session = session, inputId = "circBorderCol5", value = "#000000")
    colourpicker::updateColourInput(session = session, inputId = "labCol5", value = "#000000")
    updateSelectInput(session = session, inputId = "fill5", selected = "d\\")
    updateNumericInput(session = session, inputId = "nFills5", value = "3")
    colourpicker::updateColourInput(session = session, inputId = "fillCol5", value = "#000000")
    
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
#     colourpicker::updateColourInput(session = session, inputId = "labCol1", value = "#000000")
#     colourpicker::updateColourInput(session = session, inputId = "labCol2", value = "#000000")
#     colourpicker::updateColourInput(session = session, inputId = "labCol3", value = "#000000")
#     colourpicker::updateColourInput(session = session, inputId = "labCol4", value = "#000000")
#     colourpicker::updateColourInput(session = session, inputId = "labCol5", value = "#000000")
    
  })
  
  observeEvent(input$style_round, {
    if (input$innRatio != 0) {
      updateNumericInput(session, "innRatio", value = 0)
    } else {
      updateNumericInput(session, "innRatio", value = 0.4)  
    }
  })
  
  observeEvent(input$style_pattern, {
    
    if (any(c(input$fill1, input$fill2, input$fill3, input$fill4, input$fill5) %in% c("v", "h", "d/", "d\\"))) {
      updateSelectInput(session = session, inputId = "fill1", selected = "n")
      updateSelectInput(session = session, inputId = "fill2", selected = "n")
      updateSelectInput(session = session, inputId = "fill3", selected = "n")
      updateSelectInput(session = session, inputId = "fill4", selected = "n")
      updateSelectInput(session = session, inputId = "fill5", selected = "n")
    } else {
      updateSelectInput(session = session, inputId = "fill1", selected = "n")
      updateSelectInput(session = session, inputId = "fill2", selected = "v")
      updateSelectInput(session = session, inputId = "fill3", selected = "h")
      updateSelectInput(session = session, inputId = "fill4", selected = "d/")
      updateSelectInput(session = session, inputId = "fill5", selected = "d\\") 
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
  
  observeEvent(input$style_netLeg0, {
    updateSelectInput(session, "showLeg", selected = "No")
    updateNumericInput(session, "netWidth", value = 0)
    updateNumericInput(session, "netProp", value = 2.8)
    updateNumericInput(session, "netPadL", value = 0.1)
    updateNumericInput(session, "netPadR", value = 0.1)
    updateNumericInput(session, "netPadTop", value = 0)
  })
  
  observeEvent(input$style_netLeg1, {
    updateSelectInput(session, "showLeg", selected = "Yes")
    updateNumericInput(session, "netWidth", value = 0)
    updateNumericInput(session, "netProp", value = 3.8)
    updateNumericInput(session, "netPadL", value = 0.1)
    updateNumericInput(session, "netPadR", value = 0.1)
    updateNumericInput(session, "legX", value = -0.1)
    updateNumericInput(session, "legY", value = 0.88)
    updateNumericInput(session, "netPadTop", value = 2.2)
  })
  
  observeEvent(input$style_netLeg2, {
    updateSelectInput(session, "showLeg", selected = "Yes")
    updateNumericInput(session, "netWidth", value = 240)
    updateNumericInput(session, "netProp", value = 2.8)
    updateNumericInput(session, "netPadL", value = 0.1)
    updateNumericInput(session, "netPadR", value = 0.1)
    updateNumericInput(session, "legX", value = 0.8)
    updateNumericInput(session, "legY", value = 0.9)
    updateNumericInput(session, "netPadTop", value = 0)
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
  
#   observeEvent(input$bt_rdata, {
#     save.image("Image.RData")
#   })
#   
#   observeEvent(input$bt_rdataload, {
#     load("Image.RData")
#   })

  observeEvent(input$alphaHelix, {
    updateNumericInput(session, "netPerTurn", value = 3.6)
#     updateNumericInput(session, "netDiameter", value = 4.6)
#     updateNumericInput(session, "netPitch", value = 5.4)
#     updateNumericInput(session, "netTrans", value = 1.5)
  })
 
  observeEvent(input$three10Helix, {
    updateNumericInput(session, "netPerTurn", value = 3)
#     updateNumericInput(session, "netDiameter", value = 3.8)
#     updateNumericInput(session, "netPitch", value = 6)
#     updateNumericInput(session, "netTrans", value = 2)
  })
  
  observeEvent(input$piHelix, {
    updateNumericInput(session, "netPerTurn", value = 4.4)
#     updateNumericInput(session, "netDiameter", value = 5.6)
#     updateNumericInput(session, "netPitch", value = 4.8)
#     updateNumericInput(session, "netTrans", value = 1.1)
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
      paste(fnameNet, tolower(input$expFormat), sep = ".")
    },
    content = function(file) {
      dpi <- input$expDPI
      dims <- as.numeric(dpi)*c(getWidth(), getHeight())/72
      ImgFormat <- input$expFormat
      if (ImgFormat == "PNG") {
        png(file, res = dpi, width = dims[1], height = dims[2])
      } else if (ImgFormat == "TIFF") {
        tiff(file, compression = "lzw", res = dpi, width = dims[1], height = dims[2])
      } else if (ImgFormat == "PDF") {
        pdf(file, width = getWidth()/72, height = getHeight()/72)
      } else if (ImgFormat == "JPG") {
        jpeg(file, res = dpi, width = dims[1], height = dims[2], quality = 100)
      }
      par(mar=c(0,0,0,0) + 0.1)
      netPlot()
      dev.off()
    }
  )
})
