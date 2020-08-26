package clas12.mon.exclusive

import org.jlab.clas.physics.LorentzVector
import org.jlab.groot.data.H1F
import org.jlab.groot.data.H2F
import java.util.concurrent.ConcurrentHashMap


class EPPipPim_mon { //defining a relevant class
  def hists = new ConcurrentHashMap() 

  def beam = LorentzVector.withPID(11,0,0,10.6) //Uses pid to create electron bean lorentz vector
  def target = LorentzVector.withPID(2212,0,0,0) //the target is a stationary proton in the lab frame

  def hmm2 = {new H1F("$it","$it",250,-0.5,2)} //creating mass squared histogram
  def him = {new H1F("$it","$it",180,0.2,2)} //creating mass histogram

  def banknames = ['REC::Particle'] //identifies particle bank
  def processEvent(event) { //new function definition
    if(banknames.every{event.hasBank(it)}) {  //Checks the banknames list for elements which contain banks
      def (partb) = banknames.collect{event.getBank(it)} //loops through the events of the current bank

      (0..<partb.rows()).findAll{partb.getInt('pid',it)==11 && partb.getShort("status",it)<0} //collects all relevant events, storing them in in a list (electron, proton, pi+, and pi-)
        .collectMany{iele-> //" "
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==2212}.collect{ipro->[iele,ipro]} 
        }.collectMany{iele,ipro->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==211}.collect{ipip->[iele,ipro,ipip]}
        }.collectMany{iele,ipro,ipip->
          (0..<partb.rows()).findAll{partb.getInt('pid',it)==-211}.collect{ipim->[iele,ipro,ipip,ipim]}
        }.each{iele,ipro,ipip,ipim->
          def ele = LorentzVector.withPID(11,*['px','py','pz'].collect{partb.getFloat(it,iele)}) //reconstructing a post-scattering lorentz vector for each particle
          def pro = LorentzVector.withPID(2212,*['px','py','pz'].collect{partb.getFloat(it,ipro)})
          def pip = LorentzVector.withPID(211,*['px','py','pz'].collect{partb.getFloat(it,ipip)})
          def pim = LorentzVector.withPID(-211,*['px','py','pz'].collect{partb.getFloat(it,ipim)})
          def mm2pip = beam+target-ele-pro-pim //Calculating (lorentz vectors which correspond to) squared masses
          def mm2pim = beam+target-ele-pro-pip
          def mm2pro = beam+target-ele-pip-pim
          def impropip = pro+pip 
          def impropim = pro+pim
          def impippim = pip+pim

          def prodet = (partb.getShort('status',ipro)/1000).toInteger()==2 ? 'FD':'CD' 
          def pipdet = (partb.getShort('status',ipip)/1000).toInteger()==2 ? 'FD':'CD'
          def pimdet = (partb.getShort('status',ipim)/1000).toInteger()==2 ? 'FD':'CD'

          hists.computeIfAbsent("mm2pip_$pimdet",hmm2).fill(mm2pip.mass2()) //filling histograms with mass and mass squared values
          hists.computeIfAbsent("mm2pim_$pipdet",hmm2).fill(mm2pim.mass2())
          hists.computeIfAbsent("mm2pro_$prodet",hmm2).fill(mm2pro.mass2())
          hists.computeIfAbsent("impropip",him).fill(impropip.mass())
          hists.computeIfAbsent("impropim",him).fill(impropim.mass())
          hists.computeIfAbsent("impippim",him).fill(impippim.mass())

          [
            ['m2pip.lt.0.5/m2pim.lt.0.5', mm2pim.mass2().abs()<0.5 && mm2pip.mass2().abs()<0.5], //filling histograms with a cut
          ].findAll{it[1]}.each{name,cut->
            hists.computeIfAbsent("$name/mm2pip_$pimdet",hmm2).fill(mm2pip.mass2())
            hists.computeIfAbsent("$name/mm2pim_$pipdet",hmm2).fill(mm2pim.mass2())
            hists.computeIfAbsent("$name/mm2pro_$prodet",hmm2).fill(mm2pro.mass2())
            hists.computeIfAbsent("$name/impropip",him).fill(impropip.mass())
            hists.computeIfAbsent("$name/impropim",him).fill(impropim.mass())
            hists.computeIfAbsent("$name/impippim",him).fill(impippim.mass())
          }
        }
    }
  }
}
