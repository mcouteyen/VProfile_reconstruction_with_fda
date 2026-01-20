
data_indexes <- unique(sparse_data$INDEX)

ui <- page_sidebar(
  sidebar = sidebar(selectInput(
      inputId = "fda_index",
      label = "Index:",
      choices = data_indexes,
      selected = data_indexes[1]
    ),
    shinyWidgets::sliderTextInput(
      "fda_slider_lambda",
      "Set penalty parameter lambda",
      choices=c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2),
      selected=1e-4, 
      grid = T),
      # select K
      sliderInput(
        inputId = "fda_slider_K",
        label = "Set the number of basis functions K",
        min = 30,
        max = 250,
        value = 200
      )
    ),
    withSpinner(plotOutput("fda_plot", height = "650px"))
)

server <- function(input, output, session) {
  values <- reactiveValues()
  
  wobj_to_df <- function(W_obj){
    basis <- W_obj$basis
    n_curves <- ncol(W_obj$coefs)
    args <- seq(basis$rangeval[1], basis$rangeval[2], .2)
    df_all <- NULL
    for(i in 1:n_curves){
      profile <- fd(basisobj=basis)
      profile$coefs <- W_obj$coefs[,i]
      df_profile <- data.frame(N = i, 
                               x = args, 
                               yhat = predict(profile, newdata=args))
      df_all <- rbind(df_all, df_profile)
    }
    df_all$yhat <- df_all$yhat
    return(df_all)
  }
  
  compute_W_obj <- function(data, lambda, K, min_arg = NA, max_arg = NA){
    data <- data[!is.na(data[, my_var]),]
    #compute spline
    arg_prof <- unique(data[, "DEPTH"])
    if(is.na(min_arg)){
      min_arg <- min(arg_prof)
    }
    if(is.na(max_arg)){
      max_arg <- max(arg_prof)
    }
    n_points <- length(arg_prof)
    dates <- unique(data[, "INDEX"])
    n_curves <- length(dates)
    n_basis <- K
    
    basis <- create.bspline.basis(c(min_arg, max_arg), nbasis=n_basis, norder=4)
    
    coefs <- matrix(nrow = basis$nbasis)
    for(index in 1:n_curves){
      date_data <- data[which(data[, "INDEX"] == dates[index]),]
      date_argvals <- date_data[, "DEPTH"]
      if(min(date_argvals) <= (min_arg + 1) & max(date_argvals) >= (max_arg - 1)){
        date_y <- date_data[,my_var]
        date_smooth <- smooth.basisPar(argvals = date_argvals, y = date_y, fdobj = basis, lambda = lambda)
        coefs <- cbind(coefs, date_smooth$fd$coefs)
      }
    }
    
    W_obj <- fd(coefs, basis)
    return(W_obj)
  }
  
  plot_spline_lambda <- function(data, lambda, K){
    min_arg <- min(data[,"DEPTH"])
    max_arg <- max(data[,"DEPTH"])
    
    Wobj <- compute_W_obj(data, lambda, K, min_arg, max_arg)
    df_all <- wobj_to_df(Wobj)
    
    res <- ggplot()  +
      geom_path(data = df_all, aes(yhat, x, group=N), size=.7) +
      geom_point(data = data, aes(.data[[my_var]], DEPTH), shape=1, size=5) +
      scale_y_reverse() +
      scale_x_continuous(expand = c(0.3, 0.3)) + 
      xlab(my_var) + 
      ylab("DEPTH") + 
      theme_light() + 
      theme(legend.text = element_text(size = 24),
            legend.title = element_text(size = 24),
            axis.text.x =  element_text(size = 24),
            axis.text.y =  element_text(size = 24),
            axis.title = element_text(size=24))
    return(res)
  }
  
  output$fda_plot <- renderPlot({
    sparse_data <- sparse_data %>% filter(INDEX == input$fda_index)
    
    my_plot <- plot_spline_lambda(sparse_data, input$fda_slider_lambda, input$fda_slider_K)
    print(my_plot)
  })
}





