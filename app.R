library(shiny)
library(MASS)
library(mnormt)
library(condMVNorm)
library(mvtnorm)
library(HH)
library(plotly)
library(fMultivar)
library(tidyverse)
library(plot3D)
library(DT)
library(finalfit)
library(stargazer)
library(shinyWidgets)
library(ggpmisc)

source("global.R")

ui <- tabsetPanel(
  tabPanel("Distribution Bivariée",
  fluidPage(
  titlePanel("Distribution Normale Bivariée"),
  withMathJax(),
  tags$div(HTML("<script type='text/x-mathjax-config' >
                MathJax.Hub.Config({
                tex2jax: {inlineMath: [['$','$'], ['\\(','\\)']]}
                });
                </script >
                ")),
  sidebarLayout(
    sidebarPanel(
      fluidRow(
        column(4,
               sliderInput("corr",
        "Corrélation $(\\rho_{xy})$",
        min = -0.999,
        max = 0.999,
        step = 0.05,
        value = 0.7,
        animate = list(loop = TRUE, interval = 300)
      )),
      column(4,
             sliderInput("ltheta",
                         "Illumination $(long.)$",
                         min = -180,
                         max = 180,
                         step = 10,
                         value = 0,
                         animate = list(loop = TRUE, interval = 200)
             )
      ),
      column(4,
             sliderInput("lphi",
                         "Illumination $(lat.)$",
                         min = -180,
                         max = 180,
                         step = 10,
                         value = 0,
                         animate = list(loop = TRUE, interval = 200)
             ))),
      fluidRow(
        column(
          4,
          sliderInput("theta",
            "Longitude",
            min = 0,
            max = 180,
            step = 10,
            value = 0,
            animate = list(loop = TRUE, interval = 500)
          )
        ),
        column(
          4,
          sliderInput("phi",
            "Latitude",
            min = 0,
            max = 90,
            step = 10,
            value = 20,
            animate = list(loop = TRUE, interval = 300)
          )
        ),
        column(
          4,
          sliderInput("d",
                      "Perspective",
                      min = 0.25,
                      max = 4,
                      step = 0.05,
                      value = 1,
                      animate = list(loop = TRUE, interval = 300)
          )
        )
        
      ),
      fluidRow(
        column(
          4,
          sliderInput("expand",
            "Facteur d'expansion",
            min = 0.25,
            max = 2,
            step = 0.05,
            value = 0.45
          )
        ),
        column(
          4,
          sliderInput("shade",
            "Ombrage",
            min = 0,
            max = 1,
            step = 0.05,
            value = 0.25,
            animate = list(loop = TRUE, interval = 200)
          )
        ),
        column(
          4,
          sliderInput("r",
            "Distance",
            min = 0,
            max = 10,
            step = 1,
            value = 1,
            animate = list(loop = TRUE, interval = 200)
          )
        )
      ),
      
      fluidRow(
        column(
          4,
          checkboxInput("box",
            "Tracer la boîte",
            value = TRUE
          )
        ),
        column(
          4,
          conditionalPanel(
            condition = "input.box == TRUE",
            radioButtons("tick",
              "Axes",
              choices = c("détaillés" = "detailed", "simples" = "simple"),
              selected = "simple",
              inline = FALSE
            )
          )
        ),
        column(
          4,
          colourpicker::colourInput("col",
            "Couleur",
            "lightblue",
            showColour = "background",
            allowTransparent = TRUE)
        )
      ), p(),
      wellPanel(style = "background: lightyellow",
                helpText(        
                  h4("Distribution ", strong("bivariée:")), 
                  h5(
                   "$$  f_{(x_1, x_2)}=\\frac{1}{2 \\pi \\sigma_1 \\sigma_2 \\sqrt{1-\\rho^2}} e^{\\Big\\{-\\frac{1}{2(1-\\rho^2)} \\left [ z_1^2 -2 \\rho  z_1 z_2 +z_2^2  \\right ]\\Big\\}}$$",
                   "$$z_1=\\frac{x_1-\\mu_1}{\\sigma_1} \\; et \\; z_1=\\frac{x_2-\\mu_2}{\\sigma_2}  $$"),
                  h4("Distribution ", strong("multivariée:")), 
                  h5(    
                   "$$f(\\boldsymbol{x_1},\\boldsymbol{x_2} ) = \\frac{1}{2 \\pi |\\Sigma|^{1/2}} e^{\\Big\\{ -\\frac{1}{2} (\\boldsymbol{x} - \\boldsymbol{\\mu})^{\\top} \\Sigma^{-1} (\\boldsymbol{x} - \\boldsymbol{\\mu}) \\Big\\}} $$"
                  )
                )
      )
    ),
    mainPanel(
      fluidRow(
        column(12, align = "center",
               plotOutput("distPlot", width = "100%") 
               )
        ),
      fluidRow(
        column(12, align = "center",
               plotlyOutput("contour", width = "100%")
               )
        ),

    )
  )
)),
tabPanel("Distribution Conditionnelle",
         sidebarLayout(
         sidebarPanel(
           sliderInput("rhoxy",
                                "$\\rho_{xy}$",
                                min = -1,
                                max = 1,
                                step = 0.01,
                                value = 0,
                                animate = list(loop = TRUE, interval = 200)

         ),
           
           sliderInput("xval",
                       "X = ",
                       min = -3,
                       max = 3,
                       step = 0.01,
                       value = 0,
                       animate = list(loop = TRUE, interval = 200)

         ),
         wellPanel(style = "background: lightyellow",
         h4("Aire sous la distribution conditionnelle:"),
         radioButtons("tail",
                      "Zone étudiée:",
                      choices = list(
                        "Supérieure" = "upper",
                        "Inférieure" = "lower",
                        "Entre 2 points" = "two"
                      ),
                      selected = "lower",
                      inline = TRUE
         ),
         uiOutput("vxctrl")
         )
         ),

         mainPanel(
           plotOutput("conddist", width = "100%"), 
           plotOutput("conddist2", width = "100%")
           )
         )
         ),

tabPanel("Simulation",
         sidebarLayout(
           sidebarPanel(
             
                        fluidRow(
                          column(6,
                                 sliderInput("xmean",
                                    "$\\mu_x = $",
                                    min = 0,
                                    max = 100,
                                    step = 1,
                                    value = 0,
                                    animate = list(loop = TRUE, interval = 200))
                                 ),
                          column(6,
                                 sliderInput("xsd",
                                    "$\\sigma_x = $",
                                    min = 1,
                                    max = 20,
                                    step = 0.1,
                                    value = 1,
                                    animate = list(loop = TRUE, interval = 200))
                                 )
                          ),
                        
                        fluidRow(
                          column(6,
                                 sliderInput("ymean",
                                    "$\\mu_y = $",
                                    min = 0,
                                    max = 100,
                                    step = 1,
                                    value = 0,
                                    animate = list(loop = TRUE, interval = 200))
                                 ),
                          column(6,
                                 sliderInput("ysd",
                                    "$\\sigma_y = $",
                                    min = 1,
                                    max = 20,
                                    step = 0.1,
                                    value = 1,
                                   animate = list(loop = TRUE, interval = 200))
                                 )
                          ),
                        fluidRow(
                          column(6,
                        sliderInput("rhoxySmpl",
                                    "$\\rho_{xy}$",
                                    min = -1,
                                    max = 1,
                                    step = 0.01,
                                    value = 0,
                                    animate = list(loop = TRUE, interval = 200))
                                    
                        ),
                        column(6,
                               sliderInput("n",
                                           "n",
                                           min = 10,
                                           max = 1000,
                                           value = 200,
                                           step = 10))
                        ),
                        fluidRow(
                          column(6, offset = 3,
                               actionBttn("sample", "Échantillonner",
                                          color = "success",
                                          style = "stretch",
                                          icon = icon("gear"),
                                          block = TRUE))
           ),
           wellPanel(style = "background: lightyellow",
                     h4("Aire sous la distribution conditionnelle:"),
                     radioButtons("tail2",
                                  "Zone étudiée:",
                                  choices = list(
                                    "Supérieure" = "upper",
                                    "Inférieure" = "lower",
                                    "Entre 2 points" = "two"
                                  ),
                                  selected = "lower",
                                  inline = TRUE
                     ),
                     uiOutput("vxctrl2")
           )
           ),
           mainPanel(
            plotOutput("rawPlot", width = "100%"),
            fluidRow(
              column(5,
                     uiOutput("xval")),
              column(7,
                     uiOutput("reg"))),
            fluidRow(
              column(4,
                     wellPanel(
                     dataTableOutput("summary"))),
              column(8,
                     plotOutput("conddist3", width = "100%"))
            ),
           )
           
           )
         )
)


server <- function(input, output) {
  
  set.seed(ceiling(runif(1, 100, 9999)))
  
  r  <- reactive(mvrnorm(input$n, c(0, 0), matrix(c(1, input$rhoxySmpl, input$rhoxySmpl, 1), 2))) %>%
    bindEvent(input$sample)

  
  output$distPlot <- renderPlot(
    {
      x <- seq(-3, 3, 0.1)
      y <- seq(-3, 3, 0.1)
      mu <- c(0, 0)
      sigma <- matrix(c(1, input$corr, input$corr, 1), nrow = 2)
      f <- function(x, y) dmnorm(cbind(x, y), mu, sigma)
      z <- outer(x, y, f)

      persp(x, y, z,
        theta = input$theta,
        phi = input$phi,
        ltheta = input$ltheta,
        lphi = input$lphi,
        expand = input$expand,
        d = input$d,
        ticktype = input$tick,
        col = input$col,
        box = input$box,
        shade = input$shade,
        r = input$r,
        main = paste("Distribution bivariée \n \u03C1 = ", input$corr),
        cex.main = 2,
        cex.lab = 1.5
      )
    },
    height = 400,
    width = 1000
  )
  
  output$contour <- renderPlotly({
    
    mu1 <- mu2 <- 0
    sigma11 <- sigma22 <- 1
    sigma12 <- input$corr 
    
    x1 <- seq(mu1-3, mu1+3, length= 100)
    x2 <- seq(mu2-3, mu2+3, length= 100)
    z <- function(x1,x2){ z <- exp(-(sigma22*(x1-mu1)^2+sigma11*(x2-mu2)^2-2*sigma12*(x1-mu1)*(x2-mu2))/(2*(sigma11*sigma22-sigma12^2)))/(2*pi*sqrt(sigma11*sigma22-sigma12^2)) }
    f <- t(outer(x1,x2,z))
    
    plot_ly(x=x1,y=x2,z=f,type = "contour", width = 1000, height = 400) %>% 
      layout(xaxis=list(title="x"),yaxis=list(title="y"),
             title = "Distribution Normale Bivariée: Tracé de Contours")
  })
  
  output$conddist <- renderPlot({
    Sigma <- matrix(c(1,input$rhoxy, input$rhoxy,1), nrow=2)
#    bxy <- input$rhoxy * input$ysd / input$xsd
#    icept <- input$ymean - bxy * input$xmean
    x <- seq(-3,3,0.01)
    maintxt <- paste("Contour de la distribution normale bivariée (\u03C1 = ", input$rhoxy, ")")
    contour(x,x,outer(x,x,function(x,y){dmvnorm(cbind(x,y),sigma=Sigma)}),
            drawlabels = FALSE,
            nlevels = 15,
            main = maintxt,
            cex.main = 2)
    abline(v=input$xval, lwd=2, lty=2, col = "darkgreen")
    
    text(x = 3,
         y = ifelse(sign(input$rhoxy >= 0), -2, 2.5),
         labels=paste("x =  ", 
                      input$xval), 
         col = "red", 
         pos = 2, 
         cex=1.5)
    text(x = 3,
         y = ifelse(sign(input$rhoxy >= 0), -2.5, 2),
         labels=paste("y' = ", round(input$rhoxy * input$xval, 2)), 
         col = "red", 
         pos = 2, 
         cex=1.5)
    y <- dnorm(x, 
               mean =  0, 
               sd = sqrt(1 - input$rhoxy^2))
    lines(y-abs(min(x)),
          x + input$rhoxy * input$xval,
          lty=1,
          lwd=2, 
          col = "green")
    abline(0, 
           input$rhoxy, 
           col="red", 
           lwd=2, 
           lty=1)
    abline(h=input$rhoxy * input$xval, 
           lwd=2, 
           lty=2, 
           col="blue")
    points(input$xval, 
           input$rhoxy * input$xval, 
           pch=19, 
           cex=2, 
           col="red")
    
  },
  height = 400,
  width = 800)
  
  output$conddist2 <- renderPlot({
    
    condmean <- round(input$rhoxy * input$xval, 3)
    se.est <- round(sqrt(1 - input$rhoxy^2), 3)
    
    if (input$tail == "upper") {
      pgnorm(input$rng, 
             mean = condmean, 
             sd = se.est, 
             tail = "upper")
    } else if (input$tail == "lower") {
      pgnorm(input$rng, 
             mean = condmean, 
             sd = se.est, 
             tail = "lower")
    } else {
        pgnorm(c(input$rng[1], 
                 input$rng[2]),
               mean = condmean, 
               sd = se.est, 
               tail = "two")
      }
  },
  height = 400,
  width = 800)
  
  output$vxctrl <- renderUI({
    condmean <- input$corr * input$xval
    se.est <- sqrt(1 - input$rhoxy^2)
      if (input$tail == "lower" || input$tail == "upper") {
        sliderInput("rng",
                    "Limite:",
                    min = round(condmean - 3 * se.est, 2),
                    max = round(condmean  + 3 * se.est, 2),
                    value = round(condmean + se.est, 2),
                    step = round((condmean + 3 * se.est) / 100, 2),
                    width = "100%"
        )
      } else {
        sliderInput("rng",
                    "Limites:",
                    min = round(condmean - 3* se.est, 2),
                    max = round(condmean + 3* se.est , 2),
                    value = c(round(condmean - se.est, 2), round(condmean + se.est , 2)),
                    step = round((condmean + 3 * se.est) / 100, 2),
                    width = "100%"
        )
    } 
  })

  
  output$densityPlot2 <- renderPlot({
    x <- seq(-3, 3, 0.1)
    y <- seq(-3, 3, 0.1)
    mu <- c(0, 0)
    sigma <- matrix(c(1, input$corr, input$corr, 1), nrow = 2)
    f <- function(x, y) dmnorm(cbind(x, y), mu, sigma)
    z <- outer(x, y, f)
    
    d <- data.frame(x, y, z)
    names(d)
    
    dens <- kde2d(d$x, d$y)
    
    plot_ly(x = dens$x,
            y = dens$y,
            z = dens$z) %>% add_surface()
    
    
    d$dens <- dmvnorm(x = d)
    p4 <- plot_ly(d, x = ~ x, y = ~ y, z = ~ dens,
                  marker = list(color = ~ dens,
                                showscale = TRUE)) %>% 
      add_markers()
    
    p4
  })
  
  output$rawPlot <- renderPlot({

    r <- r()
    colnames(r) <- c("X", "Y")
    dat <- data.frame(X = r[, 1] * input$xsd + input$xmean, Y = r[, 2] * input$ysd + input$ymean) 
    reg <- lm(Y ~ X, data = dat)
    yp <- predict(reg, newdata = data.frame(X = c(input$xv)))
  
    ggplot(dat, aes(x = X, y = Y)) +
      geom_point(aes(), size = 2) +
      geom_density_2d_filled(alpha = 0.4, color = "black") +
      geom_density_2d(colour = "darkgrey", size = 1) +
      scale_fill_brewer(palette = "Reds", direction=1) +
      geom_smooth(method = "lm", colour = "yellow", fill = "cadetblue", size = 1) +
      geom_vline(xintercept = input$xv, 
                 linetype="dashed", 
                 color = "blue", 
                 size=1) +
      geom_point(aes(x=input$xv, 
                     y=yp),
                     colour="red",
                     size = 2) +
      stat_poly_eq(aes(label = after_stat(eq.label))) +
      stat_poly_eq(label.y = 0.9) 
    
  }, width = 800, height = 400)
  

  output$summary <- renderDataTable({
    r <- r()
    X <- r[, 1] * input$xsd + input$xmean
    Y <- r[, 2] * input$ysd + input$ymean
    dat <- data.frame(X, Y) 
    model <- lm(Y ~ X)
    
    tbl <- cbind(c("Moyenne", "Médiane", "Écart-Type", "Minimum", "Maximum"),
                 c(mean(X), median(X), sd(X), min(X), max(X)),
                 c(mean(Y), median(Y), sd(Y), min(Y), max(Y)))
    tbl <- data.frame(tbl)
    colnames(tbl) <- c("Statistique ", "Variable X", "Variable Y")
    tbl %>%
      datatable(
              rownames=FALSE,
              caption = htmltools::tags$caption(
                style = 'caption-side: top; text-align: center;',
                htmltools::em(h4("Description des variables"))),
              class = "cell-border stripe",
              options = list(
                info = FALSE,
                paging = FALSE,
                searching = FALSE,
                autoWidth = TRUE,
                columnDefs = list(list(width = '20%', targets = "_all"))
              )) %>%
                formatRound(columns=c("Variable X", "Variable Y"), digits=3)
  })
  
  output$reg <- renderUI({
    r <- r()
    X <- r[, 1] * input$xsd + input$xmean
    Y <- r[, 2] * input$ysd + input$ymean
    ysd <- sd(Y)
    xsd <- sd(X)
    rxy <- cor(X, Y)
    dat <- data.frame(X, Y) 
    model <- lm(Y ~ X, data=dat)
    yp <- predict(model, newdata = data.frame(X = c(input$xv)))
    
    helpText(
      h4("Régression:"),
      withMathJax(
      h5(
        "$$Y' = ", round(model$coefficients[2], 3), "(", input$xv, ") + ", round(model$coefficients[1], 3), " = \\color{red}{", round(yp, 3), "}$$",
        "$$s_{y.x}=\\sigma_y \\sqrt{1-r^2}=", round(ysd, 3), "\\sqrt{1 - ", round(rxy, 3), "^2}=\\color{red}{", round(ysd * sqrt(1 - rxy^2), 3), "}$$", 

      )
    )
    )

})
  
  output$xval <- renderUI({
    r <- r()
    X <- r[, 1] * input$xsd + input$xmean
    Y <- r[, 2] * input$ysd + input$ymean
    fluidRow(
      column(6, offset = 3,
    sliderInput("xv",
                "Valeur de X:",
                min = round(min(X), 2),
                max = round(max(X), 2),
                value = round(mean(X), 2)
      )
    )
    )
  })
  
  output$conddist3 <- renderPlot({
    
    r <- r()
    X <- r[, 1] * input$xsd + input$xmean
    Y <- r[, 2] * input$ysd + input$ymean
    ysd <- sd(Y)
    xsd <- sd(X)
    rxy <- cor(X, Y)
    dat <- data.frame(X, Y) 
    model <- lm(Y ~ X, data=dat)
    condmean <- predict(model, newdata = data.frame(X = c(input$xv)))
    se.est <- sd(Y) * sqrt(1 - cor(X, Y)^2)
    
    
    if (input$tail2 == "upper") {
      pgnorm(input$rng2, 
             mean = round(condmean, 2), 
             sd = round(se.est, 2), 
             tail = "upper")
    } else if (input$tail2 == "lower") {
      pgnorm(input$rng2, 
             mean = round(condmean, 2), 
             sd = round(se.est, 2), 
             tail = "lower")
    } else {
      pgnorm(c(input$rng2[1], 
               input$rng2[2]),
             mean = round(condmean, 2), 
             sd = round(se.est, 2), 
             tail = "two")
    }
  },
  height = 350,
  width = 500)
    
  
  output$vxctrl2 <- renderUI({
    r <- r()
    X <- r[, 1] * input$xsd + input$xmean
    Y <- r[, 2] * input$ysd + input$ymean
    ysd <- sd(Y)
    xsd <- sd(X)
    rxy <- cor(X, Y)
    dat <- data.frame(X, Y) 
    model <- lm(Y ~ X, data=dat)
    condmean <- predict(model, newdata = data.frame(X = c(input$xv)))
    se.est <- sd(Y) * sqrt(1 - cor(X, Y)^2)
    
    if (input$tail2 == "lower" || input$tail2 == "upper") {
      sliderInput("rng2",
                  "Limite:",
                  min = round(condmean - 3 * se.est, 2),
                  max = round(condmean  + 3 * se.est, 2),
                  value = round(condmean + se.est, 2),
                  step = round((condmean + 3 * se.est) / 100, 2),
                  width = "100%"
      )
    } else {
      sliderInput("rng2",
                  "Limites:",
                  min = round(condmean - 3* se.est, 2),
                  max = round(condmean + 3* se.est , 2),
                  value = c(round(condmean - se.est, 2), round(condmean + se.est , 2)),
                  step = round((condmean + 3 * se.est) / 100, 2),
                  width = "100%"
      )
    } 
  })

}

shinyApp(ui = ui, server = server)
