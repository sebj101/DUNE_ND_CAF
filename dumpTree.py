#!/usr/bin/env python

import sys
import os.path
import os
import ROOT
from optparse import OptionParser
from array import array

lar_active_vols = [ "LArActive", "PixelPlane", "LArCathode", "volLightUsPlane", "volLightDsPlane", "volLArLight", "volResistiveWire", "ResistiveField", "LArBot", "LArSubModule", "ArgonCubeActive" ]

def loop( events, tgeo, tout, nfiles, okruns ):

    if len(okruns) == 0:
        print "There are no runs in this TTree...skipping!"
        return

    print "Inside event loop with %d files and first run %d" % (nfiles, okruns[0])

    # updated geometry with less steel
    offset = [ 0., 305., 5. ]
    fvLo = [ -300., -100., 50. ]
    fvHi = [ 300., 100., 350. ]
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
        if ient % 100 == 0:
            print "Event %d of %d..." % (ient,N)
        events.GetEntry(ient)
        for ivtx,vertex in enumerate(event.Primaries):

            ## initialize output variables
            fileidx = ient/evt_per_file
            t_ifileNo[0] = okruns[fileidx]
            t_ievt[0] = ient%evt_per_file;
            t_vtx[0]=0.0; t_vtx[1]=0.0; t_vtx[2]=0.0;
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
            ## done

            # now ID numucc
            reaction=vertex.Reaction        

            # set the vertex location for output
            for i in range(3): 
                t_vtx[i] = vertex.Position[i] / 10. - offset[i] # cm

            # fiducial vertex cut
            fvCut = False
            for i in range(3):
                if t_vtx[i] < fvLo[i] or t_vtx[i] > fvHi[i]:
                    fvCut = True
            vtxv = ROOT.TVector3( t_vtx[0], t_vtx[1], t_vtx[2] )
            if fvCut:
                continue

            ileptraj = -1
            nfsp = 0
            # get the lepton kinematics from the edepsim file
            fsParticleIdx = {}
            for ipart,particle in enumerate(vertex.Particles):
                e = particle.Momentum[3]
                p = (particle.Momentum[0]**2 + particle.Momentum[1]**2 + particle.Momentum[2]**2)**0.5
                m = (e**2 - p**2)**0.5
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
                    for v in lar_active_vols:
                        if v in volName:
                            active = True
                            break
                    if active:
                        t_muonExitPt[0] = pt.X() / 10. - offset[0]
                        t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                        t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                        continue
                    # first hit outside active LAr -- determine exit
                    t_muonExitPt[0] = pt.X() / 10. - offset[0]
                    t_muonExitPt[1] = pt.Y() / 10. - offset[1]
                    t_muonExitPt[2] = pt.Z() / 10. - offset[2]
                    t_muonExitMom[0] = p.Momentum.x()
                    t_muonExitMom[1] = p.Momentum.y()
                    t_muonExitMom[2] = p.Momentum.z()
                    if abs(pt.X() / 10. - offset[0]) > 350.: muexit = 1 # side exit
                    elif abs(pt.Y() / 10. - offset[1]) > 150.: muexit = 2 # top/bottom exit
                    elif pt.Z() / 10. - offset[2] < 0.: muexit = 3 # upstream exit
                    elif pt.Z() / 10. - offset[2] > 500.: muexit = 4 # downstream exit
                    else:
                        print "Hit in %s at position (%1.1f, %1.1f, %1.1f) unknown exit!" % (volName, pt.X()/10.-offset[0], pt.Y()/10.-offset[1], pt.Z()/10.-offset[2])
                    exitP = p.Momentum
                    exitKE = (exitP.x()**2 + exitP.y()**2 + exitP.z()**2 + 105.3**2)**0.5 - 105.3
                    break

                endpt = leptraj.Points[-1].Position

                node = tgeo.FindNode( endpt.X(), endpt.Y(), endpt.Z() )

                t_lepDeath[0] = endpt.X()/10. - offset[0]
                t_lepDeath[1] = endpt.Y()/10. - offset[1]
                t_lepDeath[2] = endpt.Z()/10. - offset[2]

                endVolName = node.GetName()

                # dipole+HPGTPC
                if "volWorld" in endVolName or "volDetEnclosure" in endVolName: endVolIdx = 0 # outside detector components
                elif "volLArActive" in endVolName or "volPixelPlane" in endVolName: endVolIdx = 1 # active LAr
                elif "volLAr" in endVolName or "DsPlane" in endVolName or "UsPlane" in endVolName: endVolIdx = 2 # Passive component of LAr
                elif "volCylinder" in endVolName or "ArgonCube" in endVolName: endVolIdx = 2 # Passive component of LAr
                elif "volInsulation" in endVolName or "volGRE" in endVolName or "volSSMemb" in endVolName or "volReinforced" in endVolName: endVolIdx = 2
                elif "TPCChamber" in endVolName: endVolIdx = 6 # use "endcap yoke" idx for pressure vessel
                elif "TPC" in endVolName: endVolIdx = 3 # very rare active tpc stopper
                elif "ECALLeft" in endVolName or "ECALRight" in endVolName: endVolIdx = 4 # "endcap" ECALs
                elif "ECAL" in endVolName or "SB" in endVolName: endVolIdx = 5 # "barrel ECALs
                elif "Yoke" in endVolName: endVolIdx = 7
                elif "Mag" in endVolName: endVolIdx = 8

                # look for muon hits in the gas TPC
                hits = []
                for key in event.SegmentDetectors:
                    if key.first in ["TPC1", "TPC2"]:
                        hits += key.second

                tot_length = 0.0
                for hit in hits:
                    if hit.PrimaryId == ileptraj: # hit is due to the muon
                        # TG4HitSegment::TrackLength includes all delta-rays, which spiral in gas TPC and give ridiculously long tracks
                        hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                        hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                        tot_length += (hStop-hStart).Mag()

                t_muGArLen[0] = tot_length

                # muon reconstruction method
                # 1 = contained
                if endVolIdx == 1:
                    t_muonReco[0] = 1
                    if muexit != 0: print "Muon exit %d but end vol is active" % muexit
                # 2 = gas TPC match
                elif tot_length > 0.:
                    t_muonReco[0] = 2
                # 3 = ECAL stopping
                elif endVolIdx == 4 or endVolIdx == 5: 
                    t_muonReco[0] = 3
                # 4 = magnet/coil stopper
                elif endVolIdx == 7 or endVolIdx == 8: 
                    t_muonReco[0] = 4
                # PV
                elif endVolIdx == 6:
                    t_muonReco[0] = 9
                # 5 = passive Ar stopper
                elif endVolIdx == 2: 
                    t_muonReco[0] = 5
                # 6 = side-exiting
                elif muexit == 1:
                    t_muonReco[0] = 6
                # 7 = top/bottom-exiting
                elif muexit == 2:
                    t_muonReco[0] = 7
                # 8 = upstream-exiting
                elif muexit == 3:
                    t_muonReco[0] = 8

            # hadronic containment -- find hits in ArgonCube
            hits = []
            for key in event.SegmentDetectors:
                if key.first == "ArgonCube":
                    hits += key.second

            collar_energy = 0.
            total_energy = 0.
            track_length = [0. for i in range(nfsp)]
            for hit in hits:
                hStart = ROOT.TVector3( hit.Start[0]/10.-offset[0], hit.Start[1]/10.-offset[1], hit.Start[2]/10.-offset[2] )
                hStop = ROOT.TVector3( hit.Stop[0]/10.-offset[0], hit.Stop[1]/10.-offset[1], hit.Stop[2]/10.-offset[2] )
                if event.Trajectories[hit.PrimaryId].ParentId == -1:
                    track_length[fsParticleIdx[hit.PrimaryId]] += (hStop-hStart).Mag()

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
    t_ifileNo = array('i',[0])
    tout.Branch('ifileNo',t_ifileNo,'ifileNo/I')
    t_ievt = array('i',[0])
    tout.Branch('ievt',t_ievt,'ievt/I')
    t_p3lep = array('f',3*[0.0])
    tout.Branch('p3lep',t_p3lep,'p3lep[3]/F')
    t_vtx = array('f',3*[0.0])
    tout.Branch('vtx',t_vtx,'vtx[3]/F')
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
            fname = "%s/%02d/LAr.%s.%d.edepsim.root" % (args.topdir, run/1000, neutrino, run)
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
    print "got %d events in %d files = %1.1f events per file" % (events.GetEntries(), nfiles, 1.0*events.GetEntries()/nfiles)
    loop( events, tgeo, tout, nfiles, sorted(okruns) )

    fout.cd()
    tout.Write()




