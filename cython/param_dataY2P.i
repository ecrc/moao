


write,"*****************************************************";
write,"Fonctions a utiliser :";
write,"  nicoparam;";
write,"  getMatrices, \"/tmp\", \"ploucplouc\"";

write,"*****************************************************";



func param_Y2P(void)
{
  
  fic=open("atm-params.txt","w+");
  write,fic,format="%s\n", "Nlayer";
  write,fic,format="%d\n", tomo.learn.nb_layers;
  

  write,fic,format="%s\n", "r0 @ wfs lambda";
  write,fic,linesize=500,format="%g\n", r0Vis;

  write,fic,format="%s\n", "cn2 ESO units";
  write,fic,linesize=500,format="%g ", tomo.learn.r0;
  write,fic,"";
//  write,fic,linesize=500,format="%g ", cn2*100;
//  write,fic,"";

  write,fic,format="%s\n", "h in meters";
  write,fic,linesize=500,format="%g ", tomo.learn.altitude;
  write,fic,"";
  

  write,fic,format="%s\n", "l0 in meters";
  write,fic,linesize=500,format="%g ", 10^tomo.learn.l0;

close(fic);

  write,"";
  write,"";
  write,"";
  write,"";
  write,"";

fic=open("sys-params.txt","w+");
  write,fic,format="%s\n", "tel diam in meters";
  write,fic,format="%g\n", D;
  write,fic,format="%s\n", "cent obs";
 write,fic,format="%g\n", obstru;
  write,fic,format="%s\n", "N WFS";
  write,fic,format="%d\n", tomo.naso;
  write,fic,format="%s\n", "Nssp n subaps per wfs";
  write,fic,format="%d\n", tomo.wfs.nssp(1);

  write,fic,format="%s\n", "GsAlt";
  tmp = tomo.wfs.alt_lgs;
  nn = where(tmp!=0);
  if( is_array(nn) ) {
    tmp(nn) = 1./tmp(nn);
  }
  write,fic,format="%g ", tmp;
  write,fic,"";

  write,fic,format="%s\n", "type";
  write,fic,format="%d ", tomo.wfs.type;
  write,fic,"";
  

  write,fic,format="%s\n", "nlgs";
  write,fic,format="%d\n", nbLgsWfs;

  write,fic,format="%s\n", "alphaX in arcsec //the last value is for the truth sensor, overwrite by the function matcov_init_tomo_tiled";
  write,fic,format="%g ",tomo.wfs.x;
  write,fic,"";
  write,fic,format="%s\n", "alphaY in arcsec //the last value is for the truth sensor, overwrite by the function matcov_init_tomo_tiled";
  write,fic,format="%g ",tomo.wfs.y;
  write,fic,"";

  write,fic,format="%s\n", "XPup";
  write,fic,format="%g ", array(0.0, tomo.naso);
  write,fic,"";
  
  write,fic,format="%s\n", "YPup";
  write,fic,format="%g ", array(0.0, tomo.naso);
  write,fic,"";

  write,fic,format="%s\n", "thetaML";
  write,fic,format="%g ", array(0.0, tomo.naso);
  write,fic,"";

  write,fic,format="%s\n", "thetaCam";
  write,fic,format="%g ", array(0.0, tomo.naso);
  write,fic,"";

  write,fic,format="%s\n", "sensibilite";
  write,fic,format="%g ", array(1.0, tomo.naso);
  write,fic,"";

  write,fic,format="%s\n", "tracking";
  write,fic,format="%g ", array(0.0, 3);
  write,fic,"";

  write,fic,format="%s\n", "pasDphi";
  write,fic,format="%g\n", 0.0001;

  write,fic,format="%s\n", "ncpu";
  write,fic,format="%d\n", 12;

  write,fic,format="%s\n", "lgs_cst";
  write,fic,format="%g\n", 1.0;

  write,fic,format="%s\n", "noise_var";
  write,fic,format="%g\n", noiseVar;

  write,fic,format="%s\n", "spot_width arcsec";
  if( is_void(spotWidth) )
    write,fic,format="%g\n", 1.0;
  else
    write,fic,format="%g\n", spotWidth;
  
  write,fic,format="%s\n", "lgs_alt (m)";
  write,fic,format="%g\n", lgsH;
  
  write,fic,format="%s\n", "lgs_depth (m)";
  write,fic,format="%g\n", lgsdH;
  close(fic)
}




func getMatrices(path, suffix)
{


  fits_write, path+"/cmm_"+suffix+".fits", overwrite=1, cmm;
  fits_write, path+"/ctm_"+suffix+".fits", overwrite=1, cpm;
  fits_write, path+"/ctt_"+suffix+".fits", overwrite=1, cpp;
  fits_write, path+"/R_"+suffix+".fits", overwrite=1, R;
  fits_write, path+"/cee_"+suffix+".fits", overwrite=1, eps;
  fits_write, path+"/cvv_"+suffix+".fits", overwrite=1, cvv;
  fits_write, path+"/D_"+suffix+".fits", overwrite=1, mia;
  fits_write, path+"/Dx_"+suffix+".fits", overwrite=1, mca;

}

