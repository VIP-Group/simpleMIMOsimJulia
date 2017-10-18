#= -----------------------------------------------------------------------------
Simple MIMO simulator in Julia (v0.1)
--------------------------------------------------------------------------------
Revision history:
Oct-18-2017  v0.1   cs: initial version for release
--------------------------------------------------------------------------------
(c) 2017 Christoph Studer and Emre Gonultas
e-mail: studer@cornell.edu and eg566@cornell.edu
------------------------------------------------------------------------------=#

## we use the simPkg to define types as well as detector functions

module simPkg

## First define types

type Parameter # simulation parameters
  simName::String
  runID::Int
  save::Bool
  plot::Bool
  mod::String
  Q::Int
  MT::Int
  MR::Int
  trials::Int
  SNRdB_list::Array
  detector::Array
  symbols::Array
  bits::Array
  Es::Float64
  obj::Parameter; Parameter() = (x = new(); x.obj = x;)
end

type Result # results of simulation
  SER::Array
  BER::Array
  VER::Array
  time_elapsed::Float64
  obj::Result; Result() = (x = new(); x.obj = x;)
end

## Then define MIMO detector functions

# Zero-forcing (ZF) detector
function ZF(par,H,y)
  xhat = H\y
  _,idxhat = findmin(abs.(xhat*ones(1,length(par.symbols))-ones(par.MT,1)*par.symbols).^2,2)
  bithat = par.bits[ceil.(Int,idxhat/size(idxhat,1)),:]
  return ceil.(Int,idxhat/size(idxhat,1)),reshape(bithat,par.MT,par.Q)
end

# Unbiased linear minimum mean-square error (LMMSE) detector
function LMMSE(par,H,y,N0)
  W = (ctranspose(H)*H+(N0/par.Es)*eye(par.MT))\(ctranspose(H))
  xhat = W*y
  G = real(diag(W*H))
  _,idxhat = findmin(abs.(xhat*ones(1,length(par.symbols))-G*par.symbols).^2,2)
  bithat = par.bits[ceil.(Int,idxhat/size(idxhat,1)),:]
  return ceil.(Int,idxhat/size(idxhat,1)),reshape(bithat,par.MT,par.Q)
end

# ML detector via sphere decoding
function  ML(par,H,y)
  idxML=[]
  bitML=[]
  # -- initialization
  Radius = Inf
  PA = zeros(Int,par.MT,1) # path
  ST = zeros(par.MT,length(par.symbols)) # stack

  # -- preprocessing
  Q,R = qr(H)
  y_hat = Q'*y

  # -- add root node to stack
  Level = par.MT
  ST[Level,:] = abs.(y_hat[Level]-R[Level,Level]*transpose(par.symbols)).^2
  while ( Level<=par.MT )
    minPED,idx = findmin( ST[Level,:] )
    # -- only proceed if list is not empty
    if minPED<Inf
      ST[Level,idx] = Inf # mark child as tested
      if PA[Level+1:end,1]!=[]  #don't append empty matrix
        NewPath=cat(2,idx , PA[Level+1:end,1].')
      else
        NewPath = idx
      end
      # -- search child
      if ( minPED<Radius )
        # -- valid candidate found
        if ( Level>1 )
          # -- expand this best node
          PA[Level:end,1] = NewPath
          Level = Level-1 # downstep
          DF = R[Level,Level+1:end].' * par.symbols[PA[Level+1:end,1]]
          ST[Level,:] = minPED + abs.(y_hat[Level]-R[Level,Level]*par.symbols.'-DF[1]).^2
        else
          # -- valid leaf found
          idxML = NewPath
          bitML = par.bits[idxML,:]
          # -- update radius (radius reduction)
          Radius = minPED
        end
      end
    else
      # -- no more children to be checked
      Level=Level+1
    end
  end
  return idxML,reshape(bitML,par.MT,par.Q)
end

end # module simPkg

## import packages

using PyPlot
using JLD2
import simPkg

## main simulation routine
function simpleMIMOsim(args...)

  ## either use default parameters or parameters supplied by user
  if isempty(args)

    @printf "using default simulation settings and parameters...\n"

    # set default simulation parameters
    par = simPkg.Parameter()
    par.runID = 0 # simulation ID (used to reproduce results)
    par.save = true # save results?
    par.plot = true # plot results?
    par.MR = 4 # receive antennas
    par.MT = 4 # transmit antennas (set not larger than par.MR!)
    par.mod = "16QAM" # modulation type: "BPSK","QPSK","16QAM","64QAM","8PSK"
    par.trials = 1000 # number of Monte-Carlo par.trials (transmissions)
    par.SNRdB_list = collect(10:4:42) # list of SNR [dB] values to be simulated
    par.detector = ["ZF","LMMSE","ML"] # define par.detector(s) to be simulated

    # automatically generate reasonable simualtion name
    par.simName = string("ERR_",string(par.MR),"x",string(par.MT),"_",par.mod,"_",string(par.runID))

  elseif length(args)==1

    @printf "use custom simulation settings and parameters...\n"
    par = args[1] # only argument is par structure

  else
    error("simpleMIMOsim cannot have more than one argument")
  end

  ## simulator initialization

  # initialize random seed to enable reproducibility
  srand(par.runID)

  # define Gray-mapped constellation alphabets (QAM ones are according to IEEE 802.11)
  if par.mod== "BPSK"
    par.symbols = [ -1,1 ].'
  elseif par.mod== "QPSK"
    par.symbols = [ -1-1im, -1+1im,
                    +1-1im, +1+1im ].'
  elseif par.mod== "16QAM"
    par.symbols = [ -3-3im, -3-1im, -3+3im, -3+1im,
                    -1-3im, -1-1im, -1+3im, -1+1im,
                    +3-3im, +3-1im, +3+3im, +3+1im,
                    +1-3im, +1-1im, +1+3im, +1+1im].'
  elseif par.mod== "64QAM"
    par.symbols = [ -7-7im,-7-5im,-7-1im,-7-3im,-7+7im,-7+5im,-7+1im,-7+3im,
                    -5-7im,-5-5im,-5-1im,-5-3im,-5+7im,-5+5im,-5+1im,-5+3im,
                    -1-7im,-1-5im,-1-1im,-1-3im,-1+7im,-1+5im,-1+1im,-1+3im,
                    -3-7im,-3-5im,-3-1im,-3-3im,-3+7im,-3+5im,-3+1im,-3+3im,
                    +7-7im,+7-5im,+7-1im,+7-3im,+7+7im,+7+5im,+7+1im,+7+3im,
                    +5-7im,+5-5im,+5-1im,+5-3im,+5+7im,+5+5im,+5+1im,+5+3im,
                    +1-7im,+1-5im,+1-1im,+1-3im,+1+7im,+1+5im,+1+1im,+1+3im,
                    +3-7im,+3-5im,+3-1im,+3-3im,+3+7im,+3+5im,+3+1im,+3+3im ].'
  elseif par.mod=="8PSK"
    par.symbols = [ exp(1im*2*pi/8*0), exp(1im*2*pi/8*1),
                    exp(1im*2*pi/8*7), exp(1im*2*pi/8*6),
                    exp(1im*2*pi/8*3), exp(1im*2*pi/8*2),
                    exp(1im*2*pi/8*4), exp(1im*2*pi/8*5)].'
  end
  # extract average symbol energy
  par.Es = mean(abs.(par.symbols).^2)

  # precompute bit labels
  par.Q = Int(log2(length(par.symbols))) # number of bits per symbol
  par_bits_temp = bin.(0:length(par.symbols)-1,par.Q) #q is integer
  par.bits = zeros(Int,length(par.symbols),par.Q)
  for i=1:length(par.symbols)
    par.bits[i,:] = map(parse,split.(par_bits_temp[i],""))
  end

  # initialize result arrays (par.detector x SNR)
  res = simPkg.Result()
  res.VER = zeros(length(par.detector),length(par.SNRdB_list)) # vector error rate
  res.SER = zeros(length(par.detector),length(par.SNRdB_list)) # symbol error rate
  res.BER = zeros(length(par.detector),length(par.SNRdB_list)) # bit error rate

  # track simulation time
  time_elapsed = 0
  tic_time = time() # start timer

  ## start simulation

  # generate random bit stream (antenna x bit x trial)
  bits = rand(0:1,par.MT,par.Q,par.trials)

  # trials loop
  for t=1:par.trials

    # generate transmit symbol
    idx = 1+(bits[:,:,t] * (2.^(par.Q-1:-1:0)))
    s = par.symbols[idx]

    # generate iid Gaussian channel matrix & noise vector
    n = sqrt(0.5)*(randn(par.MR,1)+1im*randn(par.MR,1))
    H = sqrt(0.5)*(randn(par.MR,par.MT)+1im*randn(par.MR,par.MT))

    # transmit over noiseless channel (will be used later)
    x = H*s

    # SNR loop
    for k=1:length(par.SNRdB_list)

      # compute noise variance (average SNR per receive antenna is: SNR=par.MT*par.Es/N0)
      N0 = par.MT*par.Es*10^(-par.SNRdB_list[k]/10)

      # transmit data over noisy channel
      y = x+sqrt(N0)*n

      # algorithm loop
      for d=1:length(par.detector)

        if par.detector[d]=="ZF" # zero-forcing detector
          idxhat,bithat = simPkg.ZF(par,H,y)
        elseif par.detector[d]=="LMMSE" # unbiased MMSE detector
          idxhat,bithat = simPkg.LMMSE(par,H,y,N0)
        elseif par.detector[d]=="ML" # ML detector via sphere decoding
          idxhat,bithat = simPkg.ML(par,H,y)
        else
          error("par.detector type not defined.")
        end

        # -- compute error statistics
        err = (idx.!=idxhat)
        res.VER[d,k] = res.VER[d,k] + any(err)
        res.SER[d,k] = res.SER[d,k] + sum(err)/par.MT
        res.BER[d,k] = res.BER[d,k] + sum(sum(bits[:,:,t].!=bithat))/(par.MT*par.Q)

      end # algorithm loop

    end # SNR loop

    # keep track of remaining simulation time
    if time()-tic_time>10 # print a message every 10 seconds
      time_elapsed += time()-tic_time # basically adds toc()
      @printf "estimated remaining simulation time: %3.1f min.\n" time_elapsed*(par.trials/t-1)/60
      tic_time = time() # new tic()
    end

  end # par.trials loop

  # normalize results
  res.VER = res.VER/par.trials
  res.SER = res.SER/par.trials
  res.BER = res.BER./par.trials
  res.time_elapsed = time_elapsed

  if par.save
    ## save final results (par and res structure)
    @save string(par.simName,".jld2") par res
  end

  if par.plot

    ## display results (generates a nice plot)
    marker_style = ["bo-","rs--","mv-.","kp:","g*-","c>--","yx:"]
    close() # closes all figures
    figure(1)
    for d=1:length(par.detector)
      if d==1
        semilogy(par.SNRdB_list,res.BER[d,:],marker_style[d],LineWidth=2,label=par.detector[d])
      else
        semilogy(par.SNRdB_list,res.BER[d,:],marker_style[d],LineWidth=2,label=par.detector[d])
      end
    end
    grid("on",which="both",ls="--")
    xlabel("average SNR per receive antenna [dB]",fontsize=13)
    ylabel("bit error rate (BER)",fontsize=13)
    axis([minimum(par.SNRdB_list),maximum(par.SNRdB_list), 1e-4, 1])
    legend(loc="upper right",fontsize=13)
    tick_params("both",labelsize=13)

    if par.save
      # save plot as pdf file
      savefig(string(par.simName,".pdf"),bbox_inches="tight")
    end

  end

end

## start simulator with default parameters and measure excecution time
@time simpleMIMOsim()
