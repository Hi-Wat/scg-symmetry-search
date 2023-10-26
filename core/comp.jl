module Compute

using LinearAlgebra,SymPy,Optim
export tensor_search,SpinspaceOperation,phys_property


struct SpinspaceOperation
    pgname::String
    pgmat::Matrix{Float64} ## in Cartesian coordinates
    trans::Vector{Float64} ## in Cartesian coordinates
    spmat::Matrix{Float64} ## in Cartesian coordinates
 end
 
struct phys_property
    name::String
    rank::Int
    spin::Bool
    inv_parity::Int
    tr_parity::Int
    symmetrize::Bool  ## supported only for rank-2 tensor
end

function rotationmat(elev::Float64,azim::Float64,angle::Float64,irr::Bool)
    x = cos(azim)*sin(elev)
    y = sin(azim)*sin(elev)
    z = cos(elev)
    m = [cos(angle)+x^2*(1-cos(angle)) x*y*(1-cos(angle))-z*sin(angle) x*z*(1-cos(angle))+y*sin(angle);
         y*x*(1-cos(angle))+z*sin(angle) cos(angle)+y^2*(1-cos(angle)) z*y*(1-cos(angle))-x*sin(angle);
         x*z*(1-cos(angle))-y*sin(angle) z*y*(1-cos(angle))+x*sin(angle) cos(angle)+z^2*(1-cos(angle))
        ]
    if irr
        m.*=-1.0
    end
    return m
end

function purespin_axisanglesearch(
    spindim::Int,spinaxis::Vector{Float64},
    threshold=1.0e-4,testnum=10)

    if spindim in [1,2]
        cost(a) = norm(spinaxis.-[cos(a[2])*sin(a[1]),sin(a[2])*sin(a[1]),cos(a[1])])
        for _ in 1:testnum
            res = optimize(cost,randn(3),LBFGS())
            param = Optim.minimizer(res)

            if Optim.minimum(res) < threshold
                elev, azim = param
                break
            end
        end
    else
        elev,azim = 0.0, 0.0
    end

    return [elev, azim] # in radian units
end

function purespinops_add!(
    spops::Vector{SpinspaceOperation},spindim::Int,spinaxis_args::Vector{Float64},
    threshold=1.0e-4)

    # spindim, angles = Sprot.axisSearch(moments)
    spops_copy = deepcopy(spops)

    if spindim==1
        for op in spops_copy
            twofold = rotationmat(spinaxis_args[1]+π/2.0,spinaxis_args[2],π/1.0,true)
            push!(spops,SpinspaceOperation(op.pgname,op.pgmat,op.trans,twofold*op.spmat))
        end
    elseif spindim==2
        for op in spops_copy
            mirror = rotationmat(spinaxis_args[1],spinaxis_args[2],π/1.0,true)
            push!(spops,SpinspaceOperation(op.pgname,op.pgmat,op.trans,mirror*op.spmat))
        end
    
    end

    return spops
end





function antiOp_transpose(resplen::Int,applen::Int)
    row = Array{Int64}([])
    col = Array{Int64}([])
    for i in 1:resplen, j in 1:applen
        push!(row,i)
        push!(col,j)
    end
    
    matlen = resplen*applen
    transpose = zeros((matlen,matlen))
    for i in 1:matlen, j in 1:matlen
        if row[i]==col[j] && col[i]==row[j]
           transpose[i,j] =  1.0
        end
    end
    return transpose
end




function symmetrize_mat(physar::Array{phys_property})
    
    mat = Vector([1.0])
    for phys in physar
        if phys.spin
            i = 0
            bufmat = [1.0]
            while i<phys.rank
                i+=1
                bufmat = kron(I(3),bufmat)
            end
            if phys.rank==2 && phys.symmetrize
                bufmat += bufmat*antiOp_transpose(3,3)
                bufmat ./=2.0
            end
            mat = kron(bufmat,mat)
        else
            i = 0
            bufmat = [1.0]
            while i<phys.rank
                i+=1
                bufmat = kron(I(3),bufmat)
            end
            if phys.rank==2 && phys.symmetrize
                bufmat += bufmat*antiOp_transpose(3,3)
                bufmat ./=2.0
            end
            mat = kron(bufmat,mat)
            
        end
    end
    
    return mat
end




function phys_operationcreate(physar::Array{phys_property},op::SpinspaceOperation,threshold=1.0e-4)
    
    mat = Vector([1.0])
    for phys in physar
        if phys.spin
            i = 0
            bufmat = [1.0]
            while i<phys.rank
                i+=1
                bufmat = kron(op.spmat,bufmat)
            end
            mat = kron(bufmat,mat)
        else
            i = 0
            bufmat = [1.0]
            while i<phys.rank
                i+=1
                bufmat = kron(op.pgmat,bufmat)
            end
            mat = kron(bufmat,mat)
            
        end
    end
    
    return mat
end




"""

interactive_bool=true :
 If true, println function to be active. default is false.

"""
function tensor_search(
    physar::Array{phys_property},operations::Vector{SpinspaceOperation},
    spindim::Int64,colaxisarg::Vector{Float64},
    appind=2,
    threshold=1.0e-4;mode="noinput",interactive_bool=false,emrep=true,collineararg_testnum=10000)
    
    if !(interactive_bool)
        oldstd = stdout
        redirect_stdout(devnull)
    end

    println("Mode is $(mode)")
    println("The array of Physical quantities is $(length(physar))")

    if mode == "eq"
        appind = length(physar)+1
        # println("Applied fields index starts from $appind")
    elseif mode == "linear_same"
        appind = length(physar)+1
        physar = hcat([physar,physar]...)
        println("Applied fields index starts from $appind")
    elseif mode == "linear"
        println("Applied fields index starts from $appind")
    else
        println("Your input cannot use...... ")
        mode = Base.prompt("Please select : 'eq' (equilibrium)/ 'linear / linear_same' (linear response) ")
    end


    respar = physar[begin:appind-1]
    appar = physar[appind:end]


    resplen = 3^(sum([resp.rank for resp in respar]))
    applen = 3^(sum([app.rank for app in appar]))
    respname = join([resp.name for resp in respar])
    appname = join([app.name for app in appar])
    println("Tensor size for response:$(resplen),  for apply:$(applen).")

    
    if mode=="eq"
        mat = zeros(Float64,(resplen*applen,resplen*applen))
    elseif mode=="linear_same"
        mat = zeros(Float64,(resplen*applen,resplen*applen))
    elseif mode=="linear"
        mat = zeros(Float64,(resplen*applen*2,resplen*applen*2)) 
    else
        @assert mode in ["eq","linear_same","linear"]
    end


    inv_parity, tr_parity, tot_tr_parity = 1, 1, 1
    
    if mode == "eq"
        for phys in physar
            tot_tr_parity *= phys.tr_parity
            if !(phys.spin)
                inv_parity *= (-1)^(phys.rank)*phys.inv_parity
                tr_parity *= phys.tr_parity
            end
        end
    elseif mode == "linear_same"
        ######### even if axial under P or T, signs are doubled due to the diagonal correlation
    elseif mode == "linear"
        for phys in respar
            tot_tr_parity *= phys.tr_parity
            if !(phys.spin)
                inv_parity *= (-1)^(phys.rank)*phys.inv_parity
                tr_parity *= phys.tr_parity
            end
        end

        for phys in appar
            tot_tr_parity *= phys.tr_parity*(-1) # to convert current conjugate to the external fields
            if !(phys.spin)
                inv_parity *= (-1)^(phys.rank)*phys.inv_parity
                tr_parity *= phys.tr_parity*(-1) # to convert current conjugate to the external fields
            end
        end
    else
        @assert mode in ["eq","linear_same","linear"]
    end





    println("Response tensor's axiality for real space under P & T operations:")
    if inv_parity==1
        inv_parity = 0
        print("P: not axial ")
    else
        print("P: axial ")
    end
    if tr_parity==1
        tr_parity = 0
        println("T: not axial ")
    else
        println("T: axial ")
    end

    


    for op in operations
        respmat = phys_operationcreate(respar,op)
        appmat = phys_operationcreate(appar,op)
        
        if abs(det(op.spmat)-1.0)<threshold
            if mode=="eq"
                mat += kron(appmat,respmat).*det(op.pgmat)^(inv_parity)
            elseif mode=="linear_same"
                mat += kron(appmat,respmat).*det(op.pgmat)^(inv_parity)
            else mode=="linear"
                mat += kron(I(2),kron(appmat,respmat)).*det(op.pgmat)^(inv_parity)

            end
        else
            if mode=="eq"
                mat += kron(appmat,respmat).*det(op.pgmat)^(inv_parity).*det(op.spmat)^(tr_parity)
            elseif mode=="linear_same"
                mat += kron(appmat,respmat)*antiOp_transpose(applen,resplen).*det(op.pgmat)^(inv_parity).*det(op.spmat)^(tr_parity)
            else mode=="linear"
                mat += kron([0.0 1.0;1.0 0.0],kron(appmat,respmat)).*det(op.pgmat)^(inv_parity).*det(op.spmat)^(tr_parity)  
            end
        end
    end
    mat ./= length(operations)
    mat[findall(x->abs(x)<threshold, mat)] .= 0.0



    if spindim==1 && norm(abs.(mat))>threshold
        matbuf = copy(mat)
        for n in 1:collineararg_testnum
            so2rotmul = SpinspaceOperation("so2",diagm([1.0,1.0,1.0]),[0.0,0.0,0.0],rotationmat(colaxisarg[1],colaxisarg[2],2.0*π/collineararg_testnum*n,false))
            respmat = phys_operationcreate(respar,so2rotmul)
            appmat = phys_operationcreate(appar,so2rotmul)
        
            if mode == "eq"
                mat += kron(appmat,respmat)*matbuf
            elseif mode == "linear_same"
                mat += kron(appmat,respmat)*matbuf
            else mode == "linear"
                mat += kron(I(2),kron(appmat,respmat))*matbuf
            end
        end
        println("=================for colinear case, SO(2) symmetry imposed=================")
        mat ./= collineararg_testnum
        mat = round.(mat,digits=3)

    end
    
    mat = round.(mat,digits=5)

    ### symmetrize the components whose "symmetrize=true"

    if mode == "eq"

        mat *= kron(symmetrize_mat(appar),symmetrize_mat(respar))
    elseif mode == "linear_same"
        mat *= kron(symmetrize_mat(appar),symmetrize_mat(respar))
    else mode == "linear"
        mat *= kron(I(2),kron(symmetrize_mat(appar),symmetrize_mat(respar)))
    end
    
        
    symbols = ["χ_{"]
    
    for (index,phys) in enumerate(physar)
        if index==appind
            symbols.*=";"
        end
        if phys.rank==0
            symbols .*= "$(phys.name)"
        else
            symbols .*= "$(phys.name)"
            for _ in 1:phys.rank
                symbols = hcat([symbols.*"$i"  for i  in 1:3]...)
            end
        end
    end
    symbols .*="}"

    if mode=="linear"
        symbolsinv = replace.(symbols,"χ"=>"κ")
        symbols = hcat(symbols...,symbolsinv...)
    end
    

    if mode == "linear" && emrep
        symbols = replace.(symbols,"χ"=>"τ^{e}")
        symbols = replace.(symbols,"κ"=>"τ^{m}")
        if tot_tr_parity==1
            println(" Forward response is χ = τ^{e}+ τ^{m}, backward response is κ = τ^{e}- τ^{m}")
            em_transmat = kron([1 1;1 -1]./sqrt(2),I(resplen*applen))
            mat = inv(em_transmat)*mat*em_transmat
            mat[findall(x->abs(x)<threshold,mat)] .=0.0
        else
            println(" Forward response is χ = τ^{e} + τ^{m}, backward response is κ = -τ^{e}+ τ^{m}")
            em_transmat = kron([1 1;-1 1]./sqrt(2),I(resplen*applen))
            mat = inv(em_transmat)*mat*em_transmat
            mat[findall(x->abs(x)<threshold,mat)] .=0.0
        end
    end

    symbolsyms = [Sym(buf) for buf in symbols] 
    println(symbolsyms)
    symbolmat = symbolsyms*mat
    tensordim = [3^phys.rank for phys in physar]

    println("")
    println("")
    println("Your input mode is $mode")
    
    if mode in ["eq","linear_same"]
        reshapedSymbolmat = reshape(symbolmat,size(zeros(tensordim...)))
#         reshapedSymbolmat = reshape(symbolmat,size(zeros(reverse(tensordim)...)))
    else
        push!(tensordim,2)
        reshapedSymbolmat = reshape(symbolmat,size(zeros(tensordim...)))
#         reshapedSymbolmat = reshape(symbolmat,size(zeros(reverse(tensordim)...)))
    end

    if !(interactive_bool)
        redirect_stdout(oldstd) # recover original stdout
    end

    return symbolmat, reshapedSymbolmat, symbolsyms
        
    
end



end
