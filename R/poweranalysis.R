#connect all libraries
create_patterns <- function(pl)
{
  #pl <- patternlayouts
  #f = 0
  patterns <- NULL
  for (f in 1:length(pl))
  {
    #f = 1

    type <- labels(pl[f])
    pd <- as.data.frame(pl[f])
    even <- 1
    pattern <- NULL
    p <- NULL
    i <- NULL
    for (i in 1:nrow(pd))
    {
      p <- array(even, c(pd[i,]))
      pattern <- append(pattern, p)

      if (even == 0)
      {
        even = 1
      }
      else if (even == 1)
      {
        even = 0
      }
    }
    patterns[[type]] <- pattern
  }
  return(patterns)
}


compare_patterns <- function(patterns, sample)
{
  #Debug
  # sample <- line$IsIntervall
  # patterns <- create_patterns(patternlayouts)
  #
  sample <-
    paste(unlist(sample, use.names = FALSE, recursive = TRUE), collapse = '')
  #Clean sample:
  sample <- str_replace_all(sample, '3', '2')
  sample <- str_replace_all(sample, '1', '0')
  sample <- str_replace_all(sample, '2', '1')
  locales <- as.data.frame(str_locate_all(sample, '1'))
  if (nrow(locales) > 2)
  {
    sample <- substr(sample, locales[1, 1], locales[nrow(locales), 1])
  }

  distances <- NULL
  i <- NULL
  patterns[['GA1']] <- array(0, nchar(sample))

  for (i in 1:length(patterns))
  {
    pattern <-
      paste(unlist(patterns[i], use.names = FALSE, recursive = TRUE),
            collapse = '')
    dist <- stringdist(sample, pattern, method = c("osa"))
    dist
    type <- labels(patterns[i])
    #distances[type] <- dist
    distances <- rbind(distances, c(type, dist))
  }
  distances

  min <- which.min(apply(distances, MARGIN = 1, min))
  predictedtype <- distances[min, 1]

  return(predictedtype)
}


predict_type <- function(types, tlength, ga3l, ga2l)
{
  idtype <- ''

  for (t in types)
  {
    if (tlength < as.integer(t[3]) && tlength > as.integer(t[2]) &&
        ga3l > as.integer(t[4]) && ga2l > as.integer(t[5]))
    {
      idtype = t[1]
      break()
    }
  }
  return(idtype)
}

engine <- function(header, data, types, t.graph, patternlayouts)
{
  header <- header
  line <- data
  #Debug
  # header = headers[[34]][[1]]
  # line = datas[[34]][[1]]

  for (i in 3:length(line))
  {
    line[, i] = as.numeric(line[, i])
  }

  line[, 11] <-
    as.double(1.1421 * line[, 4] ^ 2 - 26.69 * line[, 4] + 225.24)
  line$Power5s = as.numeric(stats::filter(
    line[, 11],
    c(0.2, 0.2, 0.2, 0.2, 0.2),
    method = c('convolution'),
    sides = 2
  ))
  line$Power10s = stats::filter(
    line[, 11],
    c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1),
    method = c('convolution'),
    sides = 2
  )
  line$Power20s = stats::filter(
    line[, 11],
    c(
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05,
      0.05
    ),
    method = c('convolution'),
    sides = 2
  )
  line$ID <- seq.int(nrow(line))

  line$Power5s[is.na(line$Power5s)] <- 0
  line$Power10s[is.na(line$Power10s)] <- 0
  line$Power20s[is.na(line$Power20s)] <- 0

  line <-
    line %>% mutate(IsIntervall = ifelse(Power20s > (mean(Power20s) + 15), 1, 0))
  nonintervallmeanpower = mean(subset(line, line$IsIntervall == 0)$Power20s)
  line <-
    line %>% mutate(IsIntervall = ifelse(Power20s > (nonintervallmeanpower + 40), 1, 0))

  intervals <- subset(line, IsIntervall == 1)
  print(paste(
    header[3][[1]][2],
    mean(intervals$HR..bpm.),
    ',',
    mean(intervals$Cadence),
    ',',
    mean(intervals$Speed..km.h.),
    ',',
    mean(intervals$Power..W.),
    ',',
    length(intervals$HR..bpm.)
  ))

  #Calculate Intervall Windows:
  tp20s <- turnpoints(line$IsIntervall)
  tps <-
    list(stop = c(as.integer(tp20s$pos), length(line$Power20s)))
  tps$start <- as.integer(c(1, as.integer(tp20s$pos)))

  #Delete empty intervalls:
  toDelete = as.vector(c())
  for (i in 1:length(tps$start))
  {
    #delete intervall when a value is NA:
    if (is.na(tps$start[i]) || is.na(tps$stop[i]))
    {
      toDelete = append(toDelete, i)
    }
    #delete intervall if shorter than 10 sec:
    else if (tps$stop[i] - tps$start[i] < 10)
    {
      toDelete = append(toDelete, i)
    }

  }
  tps$start = tps$start[-toDelete]
  tps$stop = tps$stop[-toDelete]

  #Generate intervall factor:
  tps$IntervallType <- NULL
  for (i in 1:length(tps$start))
  {
    tps$IntervallType = append(tps$IntervallType, 1)
  }

  #Determine intervall type and smooth intervalls:
  for (i in 1:length(tps$start))
  {
    #Check Non-Intervall:
    if (mean(line$Power20s[tps$start[i]:tps$stop[i]]) < (nonintervallmeanpower + 40 + 1))
    {
      line$IsIntervall[tps$start[i]:tps$stop[i]] = 1
      tps$IntervallType[i] = 1
    }
    #Check GA2:
    if (mean(line$Power20s[tps$start[i]:tps$stop[i]]) > (nonintervallmeanpower + 40) &
        mean(line$Power20s[tps$start[i]:tps$stop[i]]) < (nonintervallmeanpower + 80 + 1))
    {
      line$IsIntervall[tps$start[i]:tps$stop[i]] = 2
      tps$IntervallType[i] = 2
    }
    #Check Max Effort:
    if (mean(line$Power20s[tps$start[i]:tps$stop[i]]) > (nonintervallmeanpower + 80))
    {
      line$IsIntervall[tps$start[i]:tps$stop[i]] = 3
      tps$IntervallType[i] = 3
    }

    #smooth:
    line$IsIntervall[tps$start[i]:tps$stop[i]] = tps$IntervallType[i]
  }

  tps$IntervallType <- as.factor(tps$IntervallType)
  line$Power20s <- as.integer(line$Power20s)
  tpsrange <- range(line$Power20s)

  if (t.graph == 1)
  {
    graph <- ggplot(line)
    graph <-
      graph + geom_line(aes(y = Power20s, x = ID, colour = 'Power20s'))
    graph <-
      graph + geom_line(aes(y = HR..bpm., x = ID, colour = 'Hr'))
    graph <-
      graph + geom_line(aes(y = Speed..km.h., x = ID, colour = 'Speed'))
    graph <-
      graph + geom_rect(
        data = as.data.frame(tps),
        aes(
          xmin = start,
          xmax = stop,
          ymin = -Inf,
          ymax = +Inf,
          fill = IntervallType
        ),
        alpha = 0.2
      )
    graph
    ggsave(
      file = paste(
        gsub('[[:punct:]]', '-', header[3][[1]][2]),
        gsub(':', '-', header[4][[1]][2]),
        '_',
        header[2][[1]][2],
        '_plot.png',
        sep = ''
      ),
      plot = graph,
      dpi = 150,
      width = 16,
      height = 12,
      units = 'cm'
    )
  }
  # c(paste(
  #   mean(intervals$HR..bpm.),
  #   ',',
  #   mean(intervals$Cadence),
  #   ',',
  #   mean(intervals$Speed..km.h.),
  #   ',',
  #   mean(intervals$Power..W.),
  #   ',',
  #   length(intervals$HR..bpm.)
  # ))

  #Calculating SpeedHRFactor
  line$SpeedHRFactor <- line[, 11] / line$HR..bpm.

  #Generating Stats:
  stats.GA1 <- subset(line, line$IsIntervall == 1)
  stats.GA2 <- subset(line, line$IsIntervall == 2)
  stats.Max <- subset(line, line$IsIntervall == 3)

  c(header[3][[1]][2])

  #output <- c(as.Date(gsub('[[:punct:]], _', '-', header[3][[1]][2])), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)

  output <-
    data.frame(
      'Date' = as.Date(gsub('[[:punct:]]', '-', header[3][[1]][2])),
      'GA1_HR' = as.numeric(NA),
      'GA1_Cadence' = as.numeric(NA),
      'GA1_Speed' = as.numeric(NA),
      'GA1_Power' = as.numeric(NA),
      'GA1_Duration' = as.numeric(NA),
      'GA2_HR' = as.numeric(NA),
      'GA2_Cadence' = as.numeric(NA),
      'GA2_Speed' = as.numeric(NA),
      'GA2_Power' = as.numeric(NA),
      'GA2_Duration' = as.numeric(NA),
      'GA3_HR' = as.numeric(NA),
      'GA3_Cadence' = as.numeric(NA),
      'GA3_Speed' = as.numeric(NA),
      'GA3_Power' = as.numeric(NA),
      'GA3_Duration' = as.numeric(NA),
      TrainingType = '',
      LengthOA = as.numeric(NA),
      PredictedType = '',
      ComparedType = '',
      SpHR_OA = as.numeric(NA),
      SpHR_GA1 = as.numeric(NA),
      SpHR_GA2 = as.numeric(NA),
      SpHR_INT = as.numeric(NA)
    )

  if (exists('stats.GA1'))
  {
    output[2] <- mean(stats.GA1$HR..bpm., na.rm = TRUE)
    output[3] <- mean(stats.GA1$Cadence, na.rm = TRUE)
    output[4] <- mean(stats.GA1$Speed..km.h., na.rm = TRUE)
    output[5] <- mean(stats.GA1$Power..W., na.rm = TRUE)
    output[6] <- length(stats.GA1$HR..bpm.)
    output[22] <- mean(stats.GA1$SpeedHRFactor, na.rm = TRUE)
    output[17] <- 'GA1'
  }
  if (exists('stats.GA2'))
  {
    output[7] <- mean(stats.GA2$HR..bpm., na.rm = TRUE)
    output[8] <- mean(stats.GA2$Cadence, na.rm = TRUE)
    output[9] <- mean(stats.GA2$Speed..km.h., na.rm = TRUE)
    output[10] <- mean(stats.GA2$Power..W., na.rm = TRUE)
    output[11] <- length(stats.GA2$HR..bpm.)
    output[23] <- mean(stats.GA2$SpeedHRFactor, na.rm = TRUE)
    if (output[11] > 60)
    {
      output[17] <- 'GA2'
    }
  }
  if (exists('stats.Max'))
  {
    output[12] <- mean(stats.Max$HR..bpm., na.rm = TRUE)
    output[13] <- mean(stats.Max$Cadence, na.rm = TRUE)
    output[14] <- mean(stats.Max$Speed..km.h., na.rm = TRUE)
    output[15] <- mean(stats.Max$Power..W., na.rm = TRUE)
    output[16] <- length(stats.Max$HR..bpm.)
    output[24] <- mean(stats.Max$SpeedHRFactor, na.rm = TRUE)
    if (output[16] > 60)
    {
      output[17] <- 'Intervall'
    }
  }
  output[18] <- length(line$HR..bpm.)
  output[19] <-
    predict_type(types, output[18], output[16], output[11])
  output[20] <-
    compare_patterns(create_patterns(patternlayouts), line$IsIntervall)
  output[21] <- mean(line$SpeedHRFactor, na.rm = TRUE)


  return(output)
}

loadlibraries <- function()
{
  if (!require("pacman"))
    install.packages("pacman")
  pacman::p_load(
    plyr,
    dplyr,
    reshape,
    stringr,
    ggplot2,
    stats,
    inline,
    pastecs,
    quantmod,
    xml2,
    stringdist
  )
}

ptanalyzer <-
  function(dirin,
           dirout,
           t.graph,
           types,
           patternlayouts) {
    loadlibraries()
    tmpwd <- getwd()
    setwd(dirin)
    filenames = list.files(pattern = "*.*csv")

    datas = list()
    headers = list()

    for (i in 1:length(filenames))
    {
      headers[[i]] <-
        list(read.csv2(
          filenames[i],
          header = FALSE,
          sep = ',',
          nrows = 3,
          as.is = TRUE
        ))
      datas[[i]] <-
        list(read.csv(
          filenames[i],
          header = TRUE,
          sep = ',',
          skip = 2
        ))
    }

    setwd(dirout)

    for (i in 1:length(headers))
    {
      if (headers[[i]][[1]][2][[1]][2] == 'INDOOR_CYCLING')
      {
        if (!exists('results'))
        {
          results <-
            engine(headers[[i]][[1]], datas[[i]][[1]], types, t.graph, patternlayouts)
        }
        else
        {
          results <-
            rbind(results,
                  engine(headers[[i]][[1]], datas[[i]][[1]], types, t.graph, patternlayouts))
        }
        print(paste('File ', filenames[i], ' completed'))
      }
      else
      {
        print(paste(
          'File with unhandeled training ignored: ',
          filenames[i],
          headers[[i]][[1]][2][[1]][2]
        ))
      }
    }
    setwd(tmpwd)
    return(results)
  }


pta.facet.all.ct <- function(results)
{
  graph <- ggplot(results)
  graph <-
    graph + geom_line(aes(y = GA1_Power, x = Date, colour = 'GA1_Power20s'))
  graph <-
    graph + geom_line(aes(y = GA2_Power, x = Date, colour = 'GA2_Power20s'))
  graph <-
    graph + geom_line(aes(y = GA3_Power, x = Date, colour = 'GA3_Power20s'))
  graph <-
    graph + geom_line(aes(y = GA1_HR, x = Date, colour = 'GA1_HR'))
  graph <-
    graph + geom_line(aes(y = GA2_HR, x = Date, colour = 'GA2_HR'))
  graph <-
    graph + geom_line(aes(y = GA3_HR, x = Date, colour = 'GA3_HR'))
  graph <- graph + facet_grid(ComparedType ~ .)
  graph <- graph + ylim(100, 300)
  return(graph)
}

pta.facet.all.tt <- function(results)
{
  graph <- ggplot(results)
  graph <-
    graph + geom_point(aes(y = GA1_Power, x = Date, colour = 'GA1_Power20s')) + geom_smooth(aes(y =
                                                                                                  GA1_Power, x = Date, colour = 'GA1_Power20s'), span = 0.3)
  graph <-
    graph + geom_point(aes(y = GA2_Power, x = Date, colour = 'GA2_Power20s')) + geom_smooth(aes(y =
                                                                                                  GA2_Power, x = Date, colour = 'GA2_Power20s'), span = 0.3)
  graph <-
    graph + geom_point(aes(y = GA3_Power, x = Date, colour = 'GA3_Power20s')) + geom_smooth(aes(y =
                                                                                                  GA3_Power, x = Date, colour = 'GA3_Power20s'), span = 0.3)
  graph <-
    graph + geom_point(aes(y = GA1_HR, x = Date, colour = 'GA1_HR')) + geom_smooth(aes(y =
                                                                                         GA1_HR, x = Date, colour = 'GA1_HR'), span = 0.3)
  graph <-
    graph + geom_point(aes(y = GA2_HR, x = Date, colour = 'GA2_HR')) + geom_smooth(aes(y =
                                                                                         GA2_HR, x = Date, colour = 'GA2_HR'), span = 0.3)
  graph <-
    graph + geom_point(aes(y = GA3_HR, x = Date, colour = 'GA3_HR')) + geom_smooth(aes(y =
                                                                                         GA3_HR, x = Date, colour = 'GA3_HR'), span = 0.3)
  graph <-
    graph + facet_grid(TrainingType ~ .) + scale_x_date(date_breaks = "1 month", date_labels = "%m-%y")
  return(graph)
}

pta.facet.sphr <- function(results)
{
  graph <- ggplot(results)
  graph <-
    graph + geom_point(aes(y = SpHR_OA, x = Date, colour = 'SpHR_OA')) + geom_smooth(aes(y =
                                                                                           SpHR_OA, x = Date, colour = 'SpHR_OA'), span = 0.3)
  graph <-
    graph + geom_point(aes(y = SpHR_GA1, x = Date, colour = 'SpHR_GA1')) + geom_smooth(aes(y =
                                                                                             SpHR_GA1, x = Date, colour = 'SpHR_GA1'), span = 0.3)
  graph <-
    graph + geom_point(aes(y = SpHR_GA2, x = Date, colour = 'SpHR_GA2')) + geom_smooth(aes(y =
                                                                                             SpHR_GA2, x = Date, colour = 'SpHR_GA2'), span = 0.3)
  graph <-
    graph + geom_point(aes(y = SpHR_INT, x = Date, colour = 'SpHR_INT')) + geom_smooth(aes(y =
                                                                                             SpHR_INT, x = Date, colour = 'SpHR_INT'), span = 0.3)
  #graph <- graph + facet_grid(TrainingType ~ .)
  graph <-
    graph  + scale_x_date(date_breaks = "1 month", date_labels = "%m-%y")
  return(graph)
}

#
# graph <- ggplot(results)
# graph <- graph + geom_smooth(aes(y=GA1_Power, x=Date, colour= 'GA1_Power20s'), span = 0.3)
# graph <- graph + geom_smooth(aes(y=GA2_Power, x=Date, colour= 'GA2_Power20s'), span = 0.3)
# graph <- graph + geom_smooth(aes(y=GA3_Power, x=Date, colour= 'GA3_Power20s'), span = 0.3)
# graph <- graph + geom_smooth(aes(y=GA1_HR, x=Date, colour= 'GA1_HR'), span = 0.3)
# graph <- graph + geom_smooth(aes(y=GA2_HR, x=Date, colour= 'GA2_HR'), span = 0.3)
# graph <- graph + geom_smooth(aes(y=GA3_HR, x=Date, colour= 'GA3_HR'), span = 0.3)
# graph <- graph + facet_grid(ComparedType ~ .)
# graph <- graph + ylim(100, 300)
# graph
#
# resultsGA1 <- subset(results, TrainingType='GA1')
# resultsGA2 <- subset(results, TrainingType='GA2')
# resultsInt <- subset(results, TrainingType='Intervall')
#
# graph <- ggplot(resultsGA1)
# graph <- graph + geom_smooth(aes(y=GA1_Power, x=Date, colour= 'GA1_Power20s'), span = 0.3)
# graph <- graph + geom_smooth(aes(y=GA1_HR, x=Date, colour= 'GA1_HR'), span = 0.3)
# graph <- graph + ylim(100, 300)
# graph
#
# graph <- ggplot(resultsGA2)
# graph <- graph + geom_smooth(aes(y=GA2_Power, x=Date, colour= 'GA2_Power20s'), span = 0.3)
# graph <- graph + geom_smooth(aes(y=GA2_HR, x=Date, colour= 'GA2_HR'), span = 0.3)
# graph <- graph + ylim(100, 300)
# graph
#
# graph <- ggplot(resultsInt)
# graph <- graph + geom_line(aes(y=GA3_Power, x=Date, colour= 'Int_Power20s'))
# graph <- graph + geom_line(aes(y=GA3_HR, x=Date, colour= 'Int_HR'))
# graph <- graph + geom_smooth(aes(y=GA3_Power, x=Date, colour= 'Int_Power20s'), span = 0.2)
# graph <- graph + geom_smooth(aes(y=GA3_HR, x=Date, colour= 'Int_HR'), span = 0.3)
# graph <- graph + ylim(100, 300)
# graph
#
# graph <- ggplot(results) + geom_point(aes(x = Date, y=LengthOA))
# graph <- graph + facet_grid(TrainingType ~ .)


dirin = 'c:/Users/jensh_000/SkyDrive/Git/rScripts/Sport/data/'
dirout = 'c:/Users/jensh_000/SkyDrive/Git/rScripts/Sport/'
t.graph = 1

trainingtypes <- list(
  c('GA2+1Max', 4140, 4200, 0, 1),
  c('GA2', 5400, 5500, 0, 1),
  c('Tabata', 2340, 2420, 1, 0),
  c('7x2', 3060, 3120, 1, 0),
  c('5x4', 3330, 3390, 1, 0),
  c('5x4', 4090, 4139, 1, 0)
)

trainingpatternlayouts <- list(
  'int2x8x20' = c(
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    240,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20
  ),
  'int5x4' = c(240, 120, 240, 120, 240, 120, 240, 120, 240),
  'int5x4new' = c(240, 90, 240, 90, 240, 90, 240, 90, 240),
  'GA2' = c(1200, 300, 600, 300, 240),
  'int7x2' = c(120, 60, 120, 60, 120, 60, 120, 60, 120, 60, 120, 60, 120, 60)
)

results <- ptanalyzer(dirin, dirout, t.graph, trainingtypes, trainingpatternlayouts)

graph <- pta.facet.sphr(results)
graph


#not run

dirin = 'c:/Users/jensh_000/SkyDrive/Git/rScripts/Sport/data/'
dirout = 'c:/Users/jensh_000/SkyDrive/Git/rScripts/Sport/'
t.graph = 1

trainingtypes <- list(
  c('GA2+1Max', 4140, 4200, 0, 1),
  c('GA2', 5400, 5500, 0, 1),
  c('Tabata', 2340, 2420, 1, 0),
  c('7x2', 3060, 3120, 1, 0),
  c('5x4', 3330, 3390, 1, 0),
  c('5x4', 4090, 4139, 1, 0)
)

trainingpatternlayouts <- list(
  'int2x8x20' = c(
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    240,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20,
    10,
    20
  ),
  'int5x4' = c(240, 120, 240, 120, 240, 120, 240, 120, 240),
  'int5x4new' = c(240, 90, 240, 90, 240, 90, 240, 90, 240),
  'GA2' = c(1200, 300, 600, 300, 240),
  'int7x2' = c(120, 60, 120, 60, 120, 60, 120, 60, 120, 60, 120, 60, 120, 60)
)

results <- ptanalyzer(dirin, dirout, t.graph, trainingtypes, trainingpatternlayouts)

#not run
