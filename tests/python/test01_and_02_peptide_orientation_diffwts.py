#' Read the *_combined.dat file
#'
#' @param Workdir name and the folder name in which the _combined.dat are kept
#'
#' @return Dataframe with tilt angle, depth position, total_score
#'
#' @examples
#' read.mapping('dir', 'combined_folder')
 
def read.mapping( peptide_list, test.name ):
 
  # Create a dataframe with the energy landscapes
  main.df <- data.frame()
 
  # List out all data files in the test category
  test.path = paste( dir, test.name, sep = "" )
 
  data.dirs <- list.dirs(path = test.path)
 
 
 
  for(f in 2:length(data.dirs)) {
   
    data.files <- list.files(path = data.dirs[f], pattern = "*combined.dat")
    # Grab the pdb identifier
    #pdb.id <- strsplit(data.files[f], "_combined_ph8")[1]
    #pdb.id <- strsplit(data.files[f], "_combined_ph4")[1]
    pdb.id <- strsplit(data.files, "_combined.")[1]
    pdb.id <- unlist(pdb.id)[1]
    print(pdb.id)
    print(is.null(pdb.id))
   
    
    if(!is.null(pdb.id)){
      # Read data from landscape file
      file.name <- paste(data.dirs[f], data.files, sep = "/")
      df<- data.frame(read.table( file.name, header = T ))
     
      print(data.files)
      #df <-data.frame(read.table(file=file.name,header=T))
     
      
     
      a <- df$angle[ which(df$angle >180 ) ]
      b <- (a - 360)
      df$angle[ which(df$angle >180 ) ] <- b
     
      a <- df$azimuthal[ which(df$azimuthal >180 ) ]
      b <- (a - 360)
      df$azimuthal[ which(df$azimuthal >180 ) ] <- b
     
      df$pdb.id <- pdb.id
     
      main.df <- rbind( main.df, df )
    }
  }
  return(main.df)
}
 
#' extracts the subset of the energy landscape of a particular pdb.id
#'
#' @param df, pdb.id and the folder name in which the _combined.dat are kept
#'
#' @return Dataframe extracts energy landscape at depth=+/-15A and at azimuthal angle = 0
#'
#' @examples
#' orientation.frame(df, '1a11', 'combined_folder')
 
orientation.frame <- function( df, pdb.id, dir ) {
 
  head.df <- data.frame()
  head.df <- df[ which( df$pdb.id == pdb.id ), ]
 
  #x <- seq( head.df$zcoord[1], head.df$zcoord[nrow(head.df)], 1 )
 x <- unique( head.df$zcoord )
 
  
  y <- seq( 0,360,1 )
 
  
  z <- seq( 0,360,5 )
 
  
  zat15.df <- data.frame()
  zatm15.df <- data.frame()
  zat15.df <- head.df[ which( head.df$zcoord == 15 ), ]
  zatm15.df <- head.df[ which( head.df$zcoord == -15 ), ]
 
  min.df <- data.frame()
  previous.df <- data.frame()
  previous.df <- head.df[ which(head.df$azimuthal==0), ]
 
  dir = paste(dir,pdb.id, sep="/")
  #min.df <- read.table( file = paste(dir, paste(pdb.id,"_min_ph8.dat", sep="") ,sep="/"), header=T )
  #min.df <- read.table( file = paste(dir, paste(pdb.id,"_min_ph4.dat", sep="") ,sep="/"), header=T )
  min.df <- read.table( file = paste(dir, paste(pdb.id,"_min.dat", sep="") ,sep="/"), header=T )
 
 
  
 #_previous is at azimuthal angle = 0
  ordered.df <- previous.df[ order(previous.df$total_score), ]
  write.table( previous.df, file = paste( dir, paste(pdb.id,"_previous.dat", sep=""), sep = "/") )
 
  write.table( zat15.df, file = paste( dir, paste(pdb.id,"_at15.dat", sep=""), sep = "/") )
 
  
  write.table( zatm15.df, file = paste( dir, paste(pdb.id,"_atm15.dat", sep=""), sep = "/") )
  
}
 
#' extracts the best orientation when energy landscape
#'
#' @param df, flag=2,3 when the depths are constant; 0,1 when all variables are considered.
#'
#' @return based on the flag, it calculates the tilt, angle and depth at which there is minimum score. 
#'
#' @examples
#' compute.best.orientation(df, 0)
 
compute.best.orientation <- function( df, flag ) {
  df.out = data.frame()
  ordered.df <- df[ order(df$total_score), ]
  head <- head(ordered.df,1)
 
  if( flag == 0 || flag == 1){
   df.out = data.frame( best.x = ordered.df$angle[ which( ordered.df$total_score == head$total_score )], best.y = ordered.df$zcoord[ which( ordered.df$total_score == head$total_score )], best.z = ordered.df$azimuthal[ which( ordered.df$total_score == head$total_score )])
  } else if( flag == 2 || flag ==3 ){
   df.out = data.frame( best.x = ordered.df$angle[ which( ordered.df$total_score == head$total_score )], best.y = ordered.df$azimuthal[ which( ordered.df$total_score == head$total_score )])
  }
  
  #print(df.out)
  return(df.out)
}
 
#' plot the energy landscape
#'
#' @param dir, flag=2,3 when the depths are constant; 0,1 when all variables are considered, pdb.id
#'
#' @return based on the flag, it plots the energy landcsapes 
#'
#' @examples
#' plot.orientation.map(dir, 0, '1a11')
#'
plot.orientation.map <- function( dir, flag, pdb.id ) {
 
  
  if( flag == 0 || flag == 1){
 
    dir=paste(dir,pdb.id,sep="/")
    df_min<- data.frame(read.table( file = paste( dir, paste(pdb.id,"_min.dat", sep=""), sep = "/"), header=T ))
    df_prev <- data.frame(read.table( file = paste( dir, paste(pdb.id,"_previous.dat", sep=""), sep = "/"), header=T ))
    #total_score
    #fa_water_to_bilayer
    #f_elec_lipidlayer
    #fa_imm_elec
   
    maxlim_min <- max( df_min$total_score)
    minlim_min <- min( df_min$total_score )
    maxlim_prev <- max( df_prev$total_score )
    minlim_prev <- min( df_prev$total_score )
   
    print(maxlim_min)
    print(minlim_min)
    print(maxlim_prev)
    print(maxlim_min)
   
    
    maxlim <- maxlim_min #max(maxlim_min, maxlim_prev)
    minlim <- minlim_min #min(minlim_min, minlim_prev)
   
    if(flag==0){
      df <- df_min
     
    }else{
      df <- df_prev
    }
    #  return( min.df )
  } else if( flag == 2 || flag == 3){
    dir=paste(dir,pdb.id,sep="/")
    df_z15 <- data.frame(read.table( file = paste( dir, paste(pdb.id,"_at15.dat", sep=""), sep = "/"), header=T ))
#    df_zm15 <- data.frame( read.table( file = paste( dir, paste(pdb.id,"_atm15.dat", sep=""), sep = "/"), header=T ))
   
    maxlim_z15 <- max( df_z15$total_score )
    minlim_z15 <- min( df_z15$total_score )
#    maxlim_zm15 <- max( df_zm15$total_score )
#    minlim_zm15 <- min( df_zm15$total_score )
   
    maxlim = maxlim_z15
    minlim=minlim_z15
    #  maxlim <- max( maxlim_z15, maxlim_zm15 )
#   minlim <- min( minlim_z15, minlim_zm15 )
   
    if( flag == 2){
      df <- df_z15
    }else{
      df <- df_zm15
    }
   
  }else{
    print( "wrong choice" )
  }
 
   # Get coordinates with the best orientation
  bestorientation.df <- compute.best.orientation( df, flag )
  #print( bestorientation.df )
  df$zcoord = as.integer( df$zcoord )
if(flag==0 || flag==1)
  {p <- ggplot() +
    theme_bw() +
    geom_raster( data = df, aes( x = angle, y = zcoord, fill = total_score ) ) +
    geom_point( data = bestorientation.df, aes( x = best.x, y = best.y ), color = "white" ) +
    geom_point( data = bestorientation.df, aes( x = best.x, y = best.y ), color = "white", shape = 1, size = 0.5 ) +
    geom_point( data = experimental.data[ which( experimental.data$pdb == pdb.id ), ], aes( x = tilt, y = z_coord ), color = "red" ) +
    geom_point( data = experimental.data[ which( experimental.data$pdb == pdb.id ), ], aes( x = tilt, y = z_coord ), color = "red", shape = 2, size = 2 ) +
    scale_x_continuous( "Angle (degrees)", expand = c(0,0), breaks = c(-180, -90, 0, 90, 180) ) +
    scale_y_continuous( "zcoord (Å)", expand = c(0,0), limits = c(-60, 0) ) +
    scale_fill_viridis( "Energy\n(REU)" , limits=c( minlim, maxlim)) +
    theme( legend.position = "bottom",
           text = element_text( size = 7, color = "black"),
           axis.text = element_text( size = 7, color = "black" ),
           axis.line = element_blank(),
           axis.ticks = element_line( size = 0.35, color = "black"),
           panel.border = element_rect(),
           legend.key.size = unit(0.5, "lines"),
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
           legend.box.margin=margin(-7,-7,-2,-7))
  return(p)
}
else{
     df1 <- df[order(df$angle, df$azimuthal),]
     df=df1
    p <- ggplot( ) +
      theme_bw() +
      geom_raster( data = df, aes( x = angle, y = azimuthal, fill = total_score ) ) +
      geom_point( data = bestorientation.df, aes( x = best.x, y = best.y ), color = "white" ) +
      geom_point( data = bestorientation.df, aes( x = best.x, y = best.y ), color = "white", shape = 1, size = 0.5 ) +
      scale_x_continuous( "Normal Angle (degrees)", expand = c(0,0), breaks = c(-180, -90, 0, 90, 180) ) +
      scale_y_continuous( "Azimuthal Angle (degrees)", expand = c(0,0), breaks = c(-180, -90, 0, 90, 180) ) +
      scale_fill_viridis( "Energy\n(REU)" , limits=c( minlim, maxlim)) +
      theme( legend.position = "bottom",
             text = element_text( size = 7, color = "black"),
             axis.text = element_text( size = 7, color = "black" ),
             axis.line = element_blank(),
             axis.ticks = element_line( size = 0.35, color = "black" ),
             panel.border = element_rect(),
             legend.key.size = unit(0.5, "lines"),
             legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"),
             legend.box.margin=margin(-7,-7,-2,-7))
    return(p)
  
  }
 
 
}
#workdir <- "/path-to-workdir/"
#workdir <- "/Users/rsamant2/Dropbox/Research_gray_lab/Membrane_protein_Electrostatic_model/Benchmark_tests"
#workdir <-"/Users/rsamant2/software/Implicit_Membrane_Energy_benchmark_recovered/Implicit_Membrane_Energy_Function_Benchmark_master/data/franklin2019/tm-peptide-tilt-angle/weights_from_test7/"
#workdir <-"/Users/rsamant2/software/Implicit_Membrane_Energy_benchmark_recovered/Implicit_Membrane_Energy_Function_Benchmark_master/data/franklin2021/amphiphilic-peptide-tilt-angle/weights_from_test7_relaxed/"
workdir <-"/Users/rsamant2/software/Implicit_Membrane_Energy_benchmark_recovered/Implicit_Membrane_Energy_Function_Benchmark_master/data/mpframework_pHmode_fa_2015/tm-peptide-tilt-angle/weights_from_test7/"
 
experimental.data <- read.table( paste( workdir, "experimental_tilt.dat", sep = ""), header = T )
 
#fa_wb_1.629_felecbilayer_0.136_fimm_-1.561
#fa_wb_1.484_felecbilayer_-0.462_fimm_0.571
 
#wt_fa_water_to_bilayer = c("1.168", "1.1", "0.966", "0.672")#, 1.484]
#wt_fa_imm_elec = c("0.045", "0.009", "0.379", "0.282")#, 0.571]
#wt_f_elec_bilayer = c("0.126", "0.076", "0.016", "0.061")
wt_fa_water_to_bilayer = c("0.966")#, "0.672")#, 1.484]
wt_fa_imm_elec = c("0.379")#, "0.282")#, 0.571]
wt_f_elec_bilayer = c("0.016")#, "0.061")
 
#wt_fa_water_to_bilayer = c("1.629")#1.179","1.375","1.571")  # ,"0.90","1.1","1.3","1.5"]
#wt_f_elec_bilayer = c("0.136")#-0.107","0.286","-0.179") #,"0.05","0.07","0.10")
#wt_fa_imm_elec = c("0.126", "0.076", "0.016", "0.061")#-0.143","-0.025") #,"0.10")#,"1.0")#c("0.02","0.05","0.07","0.10")
 
main.df <-data.frame()
 
for(index1 in 1:(length(wt_fa_water_to_bilayer))) {
#  for(index2 in 1:(length(wt_f_elec_bilayer))){
#    for(index3 in 1:(length(wt_fa_imm_elec))){
          variable_folder <-paste((paste((paste((paste((paste("fa_wb_",wt_fa_water_to_bilayer[index1],sep="")),"felecbilayer_", sep="_")),wt_f_elec_bilayer[index1], sep="")),"fimm_",sep="_")),wt_fa_imm_elec[index1],sep="")
          #variable_name <- paste((paste("flh_",wt_f_elec_bilayer[j], sep="")),(paste("fde_",wt_fa_imm_elec[k],sep="")), sep="")
          #workdir <- paste(workdir1,variable_folder,sep="/")
         
          test=variable_folder
          print(test)
          new_dir = paste(workdir,variable_folder, sep="")
          print(new_dir)
          #pass the name of the 'combined_folder' which contains all pdb.id_combined.dat and pdb.id_min.dat
          df.tm.natives <- read.mapping( workdir, variable_folder )
 
 
 
          #pdb_full <- c("AA20","AA25","AA28","AA30","AA35","AA40", "1a11", "WALP23", "polyA_cappedW", "polyA_cappedY","1mp6")#, "2mag", "LK_peptide_n6", "1hu5", "1hu6", "1hu7", "1f0d", "1f0g")
          pdb_full <- c("WALP19","WALP25","WALP31","WALP35","WALP39","AA20","AA25","AA28","AA30","AA35","AA40")#, "1a11", "WALP23", "polyA_cappedW", "polyA_cappedY","1mp6","2nr1")#,"AA40"
          #pdb_full <-c('1b4v_h1','1q4g_h1','1q4g_h2','1q4g_h3','1q4g_h4','1rhz_h1','2hih_h1','2ziy_h1','3a7k_h1','3hyw_h1','3hyw_h2','3i9v_h1','3j5p_h1','3jw8_h1')
          #pdb_full <- c('1h0a_h1','1q4g_h1','1q4g_h2','1q4g_h3','1q4g_h4','1rhz_h1','2hih_h1','3hyw_h1','3hyw_h2','3i9v_h1','3j5p_h1','3jw8_h1','3tij_h1','4hhr_h1','4hhr_h2','4hhr_h3','4m5e_h1','4nwz_h1','4qnd_h1','4umw_h1','4ymk_h1','4ymk_h2','4ymk_h3','4zwn_h1','5ahv_h1','5dqq_h1','5ek8_h1','5f19_h4','5lil_h1','5mlz_h2','5uz7_h1','5w7b_h1','5w7l_h2','5w7l_h3','6an7_h1','6d26_h1','6dvy_h1','6igk_h1') 
          #pdb_full <- c("1f0d","1f0g","1hu5","1hu6","1hu7","2mag","LK_peptide_n6")#,"2mvm","2khk","5xng")
          detailed_data_full <-data.frame()
 
          numberofpdb <- length( pdb_full )
 
          for( index in 1:numberofpdb){
           
            print(index)
            pdb <- pdb_full[index]
            print(pdb)
           
            #workdir1 = paste(workdir, "adsorbed-peptide-tilt-angle" ,sep="/")
            workdir1 = paste(workdir, test, sep="")
            orientation.frame( df.tm.natives, pdb, workdir1 )
           
            
            df.tm.flag0 <- data.frame( read.table( file = paste( paste(workdir1,pdb,sep="/"), paste(pdb,"_min.dat", sep=""), sep = "/"), header=T) )
                              
            bestorientation.df <- compute.best.orientation( df.tm.flag0, 0 )
            min_tilt = bestorientation.df$best.x
            #print(min_tilt)
            min_zcoord = bestorientation.df$best.y
            #print(min_zcoord)
            min_azim = bestorientation.df$best.z
           
            p1 <- plot.orientation.map( workdir1, 0, pdb ) + ggtitle("Minimization on azimuthal angle")
           
            #df.tm.flag1 <- orientation.frame(df.tm.natives, pdb, 1)
            df.tm.flag1 <- data.frame(read.table( file = paste( paste(workdir1,pdb,sep="/"), paste(pdb,"_previous.dat", sep=""), sep = "/"), header=T ))
            bestorientation.df <- compute.best.orientation( df.tm.flag1, 1 )
            prev_tilt = bestorientation.df$best.x
            prev_zcoord = bestorientation.df$best.y
            prev_azim = bestorientation.df$best.z
           
            size_diff <- length(prev_tilt) - length(min_tilt)
            max_size <- max( length(prev_tilt), length(min_tilt) )
          
            count_loop <- 0
           
            #if( size_diff > 0 )
            #{
            #  while(count_loop < size_diff){
            
            #    min_tilt <- c( min_tilt,'NULL' )
            #    min_zcoord <- c( min_zcoord, 'NULL' )
            #    min_azim <- c( min_azim, 'NULL' )
            #    count_loop = count_loop+1
            #   }
             
            #} else if( size_diff < 0 )
            #{
            #  while(count_loop < abs(size_diff)){
               
            #    prev_tilt <- c( prev_tilt,'NULL' )
            #    prev_zcoord <- c( prev_zcoord, 'NULL' )
            #    prev_azim <- c( prev_azim, 'NULL' )
            #    count_loop = count_loop+1
           # }
             
              
           # }
           
            p2 <- plot.orientation.map( workdir1, 1, pdb) + ggtitle("At azimuthal angle = 0")
           
            exp_tilt = experimental.data[ which( experimental.data$pdb.id == pdb ), 2]
            exp_zcoord = experimental.data[ which( experimental.data$pdb.id == pdb ), 3]
           
              
            detailed_data <- data.frame( min_tilt, min_zcoord, min_azim, prev_tilt, prev_zcoord, prev_azim )
            detailed_data$pdb.id = pdb
            detailed_data$exp_tilt = exp_tilt
            detailed_data$exp_zcoord = exp_zcoord
           
              #print(detailed_data)
            detailed_data_full <- rbind(detailed_data_full, detailed_data)
            #df.tm.flag2 <- orientation.frame(df.tm.natives, pdb, 2)
            p3 <- plot.orientation.map( workdir1, 2, pdb) + ggtitle("At zcoord z = 15")
            #df.tm.flag3 <- orientation.frame(df.tm.natives, pdb, 3)
#            p4 <- plot.orientation.map( workdir1, 3, pdb) + ggtitle("At zcoord z = -15")
           
            plot_grid( p1, p2, p3, ncol = 2, nrow = 2, labels = c("a", "b", "c"), label_size = 9)
            a1.grid <- plot_grid( p1, p2,p3, ncol = 2, nrow = 2, labels = c("a", "b", "c"), label_size = 9 )
            save_plot( paste( workdir1, paste(pdb,sep="_","total_score.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 4 )
         
           # plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
           # a1.grid <- plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9 )
            #save_plot( paste( workdir, paste(pdb,sep="_","paperupdates.eps"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 4 )
           # save_plot( paste( workdir1, paste(pdb,sep="_","paperupdates.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 4 )
           
            
          }
          write.table( detailed_data_full, file = paste( new_dir, "processed_data_test1and2.txt", sep = "/") )
#    }
#  }
}
def main(args):

    all_sub_tests = ["ddG-of-insertion", "ddG-of-pH-insertion",
                     "tm-peptide-tilt-angle", "adsorbed-peptide-tilt-angle", "protein-tilt-angle"]

    # Read options from the command line
    parser = OptionParser(
        usage="usage %prog --energy_fxn franklin2019 --which_tests all")
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e', action="store",
                      help="Name of energy function weights file", )
    parser.add_option('--which_tests', '-w', action="store",
                      help="Specify tests run (comma separated list)", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check that required options have been provided
    if (not Options.energy_fxn or not Options.which_tests):
        print("Missing required options --which_tests or --energy_fxn")
        sys.exit()

    # Check test categories
    test_names = []

    # Read path configuration file
    config = read_config.read_config()

    if (Options.which_tests == "all"):
        test_names = all_sub_tests
    else:
        test_names = Options.which_tests.split(",")
        # check that all names are valid
        for name in test_names:
            if name not in all_sub_tests:
                sys.exit("No such test " + name +
                         ". Exiting! or not among the tests 1-3 and 5-6")

            else:
                
                wt_fa_water_to_bilayer = ["0.966"]#,"1.168" ]#, 1.484] 
                wt_fa_imm_elec = ["0.379"]#,"0.045"]#, 0.571]
                wt_f_elec_bilayer = ["0.016"]#,"0.126"]#, -0.462]  
        
        
                for i in range(len(wt_fa_water_to_bilayer)):
                    for j in range(len(wt_f_elec_bilayer)):
                        for k in range(len(wt_fa_imm_elec)):

                            datadir = '/home/rsamant2/scratch16-jgray21/rsamant2/' + "data/"+ \
                                Options.energy_fxn + "/" + name + "/weights_from_test7_mpenv/"
                            # datadir = config.benchmark_path + "data/" + \
                            # Options.energy_fxn + "/" + name + "/weights_from_test7/"
                            datadir = datadir + "fa_wb_" + str(wt_fa_water_to_bilayer[i]) + "_felecbilayer_" + str(
                                wt_f_elec_bilayer[i]) + "_fimm_" + str(wt_fa_imm_elec[i])
                            
                            print(datadir)
                            if (not os.path.isdir(datadir)):
                                sys.exit("No such test data available; test needs to finish before combining files")
                            os.chdir(datadir)

                            if(name == "tm-peptide-tilt-angle"):
                                # peptide_list = ['2nr1','polyA_cappedW', 'polyA_cappedY', 'AA20','AA25','AA28','AA30','AA35','AA40']
                                peptide_list = ['WALP19','WALP25','WALP31','WALP35','WALP39']
                                #peptide_list = ['AA30','AA40']
                                # peptide_list = ['AA20','AA25','AA28','AA30','AA35','AA40']#, 'WALP19','WALP25','WALP31','WALP35','WALP39']
                                # peptide_list = ['WALP25','WALP35']#['WALP19','WALP25','WALP31','WALP35','WALP39']
                            elif(name == "adsorbed-peptide-tilt-angle"):
                                peptide_list = [
                                    '1f0d', '1f0g', '1hu5', '1hu6', '1hu7', '2mag', 'LK_peptide_n6']#,'2khk','2mvm','5xng']
                            elif(name == "protein-tilt-angle"):
                                peptide_list = ['1fep', '1gzm', '1m0l', '1nqe', '1okc', '1p4t', '1qfg', '1qj8', '1qjp', '1r3j', '1yce', '2cfp', '2j8c', '2qom',
                                                '2x9k', '3aeh', '3dzm', '3pox', '3syb', '3wbn', '4afk', '4fqe', '4hyj', '4m48', '4n6h', '4rl8', '4uc2', '4x5n']
                            elif(name == "ddG-of-insertion"):
                                peptide_list = [
                                    'GL5']#'GL6', 'GL7', 'GL8', 'GWL6']#'GL5', 
                            elif(name == "ddG-of-pH-insertion"):
                                peptide_list = ['pHLIP-v1_4', 'pHLIP-v1_8', 'pHLIP-v2_4', 'pHLIP-v2_8', 'pHLIP-v3_4', 'pHLIP-v3_8', 'pHLIP-v4_4', 'pHLIP-v4_8', 'pHLIP-v5_4', 'pHLIP-v5_8', 'pHLIP-v6_4', 'pHLIP-v6_8', 'pHLIP-v7_4', 'pHLIP-v7_8', 'pHLIP-v8_4', 'pHLIP-v8_8',
                                                'pHLIP-v9_4', 'pHLIP-v9_8', 'pHLIP-v10_4', 'pHLIP-v10_8', 'pHLIP-v11_4', 'pHLIP-v11_8', 'pHLIP-v12_4', 'pHLIP-v12_8', 'pHLIP-v13_4', 'pHLIP-v13_8', 'pHLIP-v14_4', 'pHLIP-v14_8', 'pHLIP-v15_4', 'pHLIP-v15_8', 'pHLIP-v16_4', 'pHLIP-v16_8']
                            

                            

if __name__ == "__main__":
    main(sys.argv) 