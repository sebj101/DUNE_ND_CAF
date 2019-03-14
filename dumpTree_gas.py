#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array

gar_active_vols = ["GArTPC", "TPCChamber", "TPCGas", "cent_elec_shape", "cent_hc_shape", "TPC1_shape", "TPC1pad_shape", "TPC1fc_pvf_shape", "TPC1fc_kev_shape", "TPC2_shape", "TPC2pad_shape", "TPC2fc_pvf_shape", "TPC2fc_kev_shape", "TPC2fc_hc_shape"]

def loop( events, tgeo, tout, nfiles, okruns ):
 
    if len(okruns) == 0:
        print "There are no runs in this TTree...skipping!"
        return

    print "Inside event loop with %d files and first run %d" % (nfiles, okruns[0])

    # updated geometry with less steel
    offset = [ 0., 305., 5. ]
    fvLo = [ -300., -100., 50. ]
    fvHi = [ 300., 100., 450. ]
    collarLo = [ -320., -120., 30. ]
    collarHi = [ 320., 120., 470. ]

    event = ROOT.TG4Event()
    events.SetBranchAddress("Event",ROOT.AddressOf(event))
    print "Set branch address"
 
    N = events.GetEntries()
    evt_per_file = N/nfiles
    if N % nfiles:
        print "Files don't all have the same number of events!!!"
        print "AAAAAAAAAAAAAAAAAAAAAAAAHHHHHHHHHHHHHHHHhhhh that's bad stop"
        print "\n\n\n\n\n\n\n\n\n\n\n"
        print "AAAAAAAAAAAAAAAAAAAAAAAAHHHHHHHHHHHHHHHHhhhh seriously stop it"

    print "Starting loop over %d entries" % N
    for ient in range(N):
        pi_pl_count = 0
        pi_min_count = 0
       
        nue_tot = 0.0

        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)
        events.GetEntry(ient)
        for ivtx,vertex in enumerate(event.Primaries):
            node_all_prim = tgeo.FindNode(vertex.Position[0], vertex.Position[1],vertex.Position[2])
            volname = node_all_prim.GetName()
  
            fileidx = ient/evt_per_file
            t_ifileNo[0] = okruns[fileidx]
            t_ievt[0] = ient%evt_per_file;
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
                      
            t_vtx_no_off[0]=0.0; t_vtx_no_off[1]=0.0; t_vtx_no_off[2]=0.0; 
            t_p3lep[0]=0.0; t_p3lep[1]=0.0; t_p3lep[2]=0.0;
            t_lepDeath[0]=0.0; t_lepDeath[1]=0.0; t_lepDeath[2]=0.0;
            t_lepPdg[0] = 0
            t_lepKE[0] = 0.
            t_muonExitPt[0] = 0.0; t_muonExitPt[1] = 0.0; t_muonExitPt[2] = 0.0; 
            t_muonExitMom[0] = 0.0; t_muonExitMom[1] = 0.0; t_muonExitMom[2] = 0.0; 
            t_muonReco[0] = -1;
            t_muGArLen[0]=0.0;
            t_hadTot[0] = 0.
            t_hadCollar[0] = 0.
            t_nFS[0] = 0
            t_lepstart[0]=0.0; t_lepstart[1]=0.0; t_lepstart[2]=0.0;

            reaction=vertex.Reaction        

            # set the vertex location for output
            for i in range(3): 
                t_vtx[i] = vertex.Position[i] / 10. - offset[i] # cm
           
            for i in range(3):
                t_vtx_no_off[i] = vertex.Position[i] / 10.



            #for i in range(3):
            #    if t_vtx[i] < fvLo[i] or t_vtx[i] > fvHi[i]:
            #        fvCut = True
            #vtxv = ROOT.TVector3( t_vtx[0], t_vtx[1], t_vtx[2] )
            #if fvCut:
            #    continue

            ileptraj = -1
            nfsp = 0

            fsParticleIdx = {}
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                 
                p_recon = ROOT.gRandom.Gaus(p, p*(0.01))
                p_recon_pi0 = ROOT.gRandom.Gaus(p, p*(0.10))
                if particle.PDGCode != 111: 
                    t_p_recon_gen[0] = p_recon
                if particle.PDGCode == 111:
                    t_p_recon_gen[0] = p_recon_pi0
                m = (e**2 - p**2)**0.5
                t_m_recon_gen[0] = m
                t_p_true_gen[0] = p
                t_pdg_gen[0] = particle.PDGCode
        
                t_fsPdg[nfsp] = particle.PDGCode
                t_fsPx[nfsp] = particle.Momentum[0]
                t_fsPy[nfsp] = particle.Momentum[1]
                t_fsPz[nfsp] = particle.Momentum[2]
                t_fsE[nfsp] = e
                fsParticleIdx[particle.TrackId] = nfsp
                nfsp += 1

                if abs(particle.PDGCode) in [11,12,13,14]:
                    ileptraj = particle.TrackId
                    t_lepPdg[0] = particle.PDGCode
                    # set the muon momentum for output
                    for i in range(3): t_p3lep[i] = particle.Momentum[i]
                    t_lepKE[0] = e - m
                   
            assert ileptraj != -1, "There isn't a lepton??"
            t_nFS[0] = nfsp
           

            # If there is a muon, determine how to reconstruct its momentum and charge
            muexit = 0
            exitKE = 0.
            exitP = None
            endVolIdx = -1 # where does the muon die
            if abs(t_lepPdg[0]) == 13:
                leptraj = event.Trajectories[ileptraj]
                for p in leptraj.Points:
                    pt = p.Position
                    
                    node = tgeo.FindNode( pt.X(), pt.Y(), pt.Z() )
                    volName = node.GetName()
                    active = False
                    for v in gar_active_vols:
                        if v in volName:
                            active = True
                            break
                    if active:
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        continue
        
                    t_muonExitPt[0] = pt.X() / 10. - offset[0]
                    t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                    t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                    t_muonExitMom[0] = p.Momentum.x()
                    t_muonExitMom[1] = p.Momentum.y()
                    t_muonExitMom[2] = p.Momentum.z()
                    exitP = p.Momentum
                    exitKE = (exitP.x()**2 + exitP.y()**2 + exitP.z()**2 + 105.3**2)**0.5 - 105.3
                    break
                   
                  
               
                startpt = leptraj.Points[0].Position
                endpt = leptraj.Points[-1].Position
                node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

                t_lepDeath[0] = endpt.X()/10. - offset[0]
                t_lepDeath[1] = endpt.Y()/10. - offset[1]
                t_lepDeath[2] = endpt.Z()/10. - offset[2]
                
                t_lepstart[0] = startpt.X()/10. - offset[0]
                t_lepstart[1] = startpt.Y()/10. - offset[1]
                t_lepstart[2] = startpt.Z()/10. - offset[2]

                endVolName = node.GetName()


                if "volWorld" in endVolName or "volDetEnclosure" in endVolName: endVolIdx = 0 # outside detector components
                if "TPC1_PV_0" == volName or "TPC2_PV_0" == volName: 
                        endVolIdx = 1 # active GAr

                if "Endcap" in endVolName or "Yoke" in endVolName: endVolIdx = 2 # Passive component of GAr
                if "Tile" in endVolName: endVolIdx = 3

                hits = []
                
                for key in event.SegmentDetectors:
               
        
                     if key.first in ["TPC1", "TPC2"]:
                        hits += key.second

                tot_length = 0.0
                for hit in hits:
                    if hit.PrimaryId == ileptraj: # hit is due to the muon

                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        tot_length += (hStop-hStart).Mag()

                t_muGArLen[0] = tot_length

            hits = []
            for key in event.SegmentDetectors:

                if key.first in  ["TPC1", "TPC2"]:

                        hits += key.second
            collar_energy = 0.
            total_energy = 0.
            track_length = [0. for i in range(nfsp)]
            track_length_perp = [0. for i in range(nfsp)]
            for hit in hits:
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                hStartYZ = ROOT.TVector3( 0., hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStopYZ = ROOT.TVector3( 0., hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                if event.Trajectories[hit.PrimaryId].ParentId == -1:
                    track_length[fsParticleIdx[hit.PrimaryId]] += (hStop-hStart).Mag()
                    track_length_perp[fsParticleIdx[hit.PrimaryId]] += (hStopYZ-hStartYZ).Mag()
                # below is to be ignored for gas TPC
                if hit.PrimaryId != ileptraj:
                    hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                    total_energy += hit.EnergyDeposit
                    # check if hit is in collar region
                    if hStart.x() < collarLo[0] or hStart.x() > collarHi[0] or hStart.y() < collarLo[1] or hStart.y() > collarHi[1] or hStart.z() < collarLo[2] or hStart.z() > collarHi[2]:
                        collar_energy += hit.EnergyDeposit

            t_hadTot[0] = total_energy
            t_hadCollar[0] = collar_energy
            for i in range(nfsp):
                t_fsTrkLen[i] = track_length[i]
                t_fsTrkLenPerp[i] = track_length_perp[i]

            tout.Fill()


 







if __name__ == "__main__":

    ROOT.gROOT.SetBatch(1)

    parser = OptionParser()
    parser.add_option('--outfile', help='Output file name', default="out.root")
    parser.add_option('--topdir', help='Input file top directory', default="")
    parser.add_option('--first_run', type=int, help='First run number', default=0)
    parser.add_option('--last_run', type=int, help='Last run number', default=0)
    parser.add_option('--rhc', action='store_true', help='Reverse horn current', default=False)
    parser.add_option('--grid', action='store_true', help='grid mode', default=False)

    (args, dummy) = parser.parse_args()

    # make an output ntuple
    fout = ROOT.TFile( args.outfile, "RECREATE" )
    tout = ROOT.TTree( "tree","tree" )
    
   
  
 


    t_m_recon_gen = array('f',[0.0])
    tout.Branch('m_recon_gen',t_m_recon_gen,'m_recon_gen/F')
    t_pdg_gen = array('i',[0])
    tout.Branch('pdg_gen',t_pdg_gen,'pdg_gen/I') 
    t_p_true_gen = array('f',[0.0])
    tout.Branch('p_true_gen',t_p_true_gen,'p_true_gen/F') 
    t_p_recon_gen = array('f',[0.0])
    tout.Branch('p_recon_gen',t_p_recon_gen,'p_recon_gen/F')


    t_ifileNo = array('i',[0])
    tout.Branch('ifileNo',t_ifileNo,'ifileNo/I')
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')
    t_vtx_no_off = array('f',3*[0.0])
    tout.Branch('vtx_no_off',t_vtx_no_off,'vtx_no_off[3]/F')
    t_lepDeath = array('f',3*[0.0])
    tout.Branch('lepDeath',t_lepDeath,'lepDeath[3]/F')
    t_lepPdg = array('i',[0])
    tout.Branch('lepPdg',t_lepPdg,'lepPdg/I')
    t_lepKE = array('f',[0])
    tout.Branch('lepKE',t_lepKE,'lepKE/F')
    t_muonExitPt = array('f',3*[0.0])
    tout.Branch('muonExitPt',t_muonExitPt,'muonExitPt[3]/F')
    t_muonExitMom = array('f',3*[0.0])
    tout.Branch('muonExitMom',t_muonExitMom,'muonExitMom[3]/F')
    t_muonReco = array('i',[0])
    tout.Branch('muonReco',t_muonReco,'muonReco/I')
    t_muGArLen = array('f',[0])
    tout.Branch('muGArLen',t_muGArLen,'muGArLen/F')
    t_hadTot = array('f', [0.] )
    tout.Branch('hadTot', t_hadTot, 'hadTot/F' )
    t_hadCollar = array('f', [0.] )
    tout.Branch('hadCollar', t_hadCollar, 'hadCollar/F' )
    t_nFS = array('i',[0])
    tout.Branch('nFS',t_nFS,'nFS/I')
    t_fsPdg = array('i',100*[0])
    tout.Branch('fsPdg',t_fsPdg,'fsPdg[nFS]/I')
    t_fsPx = array('f',100*[0.])
    tout.Branch('fsPx',t_fsPx,'fsPx[nFS]/F')
    t_fsPy = array('f',100*[0.])
    tout.Branch('fsPy',t_fsPy,'fsPy[nFS]/F')
    t_fsPz = array('f',100*[0.])
    tout.Branch('fsPz',t_fsPz,'fsPz[nFS]/F')
    t_fsE = array('f',100*[0.])
    tout.Branch('fsE',t_fsE,'fsE[nFS]/F')
    t_fsTrkLen = array('f',100*[0.])
    tout.Branch('fsTrkLen',t_fsTrkLen,'fsTrkLen[nFS]/F')
    t_fsTrkLenPerp = array('f',100*[0.])
    tout.Branch('fsTrkLenPerp',t_fsTrkLenPerp,'fsTrkLenPerp[nFS]/F')
    t_lepstart = array('f',3*[0.0])
    tout.Branch('lepstart',t_lepstart,'lepstart[3]/F')
    
    t_nue_tot = array('f',[0.])
    tout.Branch('nue_tot',t_nue_tot,'nue_tot')
    loaded = False
    #if os.path.exists("EDepSimEvents/EDepSimEvents.so"):
    #    print "Found EDepSimEvents.so"
    #    ROOT.gSystem.Load("EDepSimEvents/EDepSimEvents.so")
    #    loaded = True

    tgeo = None

    events = ROOT.TChain( "EDepSimEvents", "main event tree" )

    neutrino = "neutrino"
    if args.rhc:
        neutrino = "antineutrino"

    print "Building TChains for runs %d-%d..." % (args.first_run, args.last_run)
    nfiles = 0
    okruns = []
    for run in range( args.first_run, args.last_run+1 ):
        if args.grid:
            fname = "%s/edep.%d.root" % (args.topdir,run)
        else:
            fname = "%s/GAr.%s.%d.edepsim.root" % (args.topdir, neutrino, run)
        print fname

        # see if it is an OK file
        if not os.access( fname, os.R_OK ):
            print "Can't access file: %s" % fname
            continue
        tf = ROOT.TFile( fname )
        if tf.TestBit(ROOT.TFile.kRecovered): # problem with file
            print "File is crap: %s" % fname
            continue
        nfiles += 1
        okruns.append( run )

        if not loaded:
            loaded = True
            tf.MakeProject("EDepSimEvents","*","RECREATE++")

        # add it to the tchain
        events.Add( fname )

        if tgeo is None: # first OK file, get geometry
            tgeo = tf.Get("EDepSimGeometry")
        tf.Close() # done with this one

    print "OK runs: ", sorted(okruns)
    #print "got %d events in %d files = %1.1f events per file" % (events.GetEntries(), nfiles, 1.0*events.GetEntries()/nfiles)
    loop( events, tgeo, tout, nfiles, sorted(okruns) )

    fout.cd()
    tout.Write()




