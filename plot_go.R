##################
## Sig of nodes ##
##################

showSigOfNodes <- function (GOdata, termsP.value, firstSigNodes = 10, reverse = TRUE, 
    sigForAll = TRUE, wantedNodes = NULL, putWN = TRUE, putCL = 0, 
    type = NULL, showEdges = TRUE, swPlot = TRUE, useFullNames = TRUE, 
    oldSigNodes = NULL, useInfo = c("none", "pval", "counts", 
        "def", "np", "all")[1], plotFunction = GOplot, .NO.CHAR = 20, sigFontSize, sigBoxSize, elseFontSize) 
{
    require("Rgraphviz") || stop("package Rgraphviz is required")
    if (!is.null(firstSigNodes)) 
        sigTerms <- sort(termsP.value)[1:firstSigNodes]
    else sigTerms <- numeric(0)
    if (putWN && !is.null(wantedNodes)) 
        baseNodes <- union(names(sigTerms), wantedNodes)
    else baseNodes <- names(sigTerms)
    if (length(baseNodes) == 0) 
        stop("No nodes were selected")
    if (putCL) {
        goDAG.r2l <- reverseArch(graph(GOdata))
        for (i in 1:putCL) {
            newNodes <- unique(unlist(adj(goDAG.r2l, baseNodes)))
            baseNodes <- union(newNodes, baseNodes)
        }
    }
    dag <- inducedGraph(graph(GOdata), baseNodes)
    if (reverse) 
        dag <- reverseArch(dag)
    termCounts <- termStat(GOdata, nodes(dag))
    if (!is.null(type)) {
        if (swPlot) 
            GOplot.counts(dag, wantedNodes = wantedNodes, nodeCounts = termCounts, 
                showEdges = showEdges)
        return(dag)
    }
    pval.info <- function(whichNodes) {
        ret.val <- format.pval(termsP.value[whichNodes], digits = 3, 
            eps = 1e-20)
        names(ret.val) <- whichNodes
        return(ret.val)
    }
    .pval = pval.info(nodes(dag))
    .def = .getTermsDefinition(whichTerms = nodes(dag), ontology(GOdata), 
        numChar = .NO.CHAR)
    .counts = apply(termCounts[, c("Significant", "Annotated")], 
        1, paste, collapse = " / ")
    nodeInfo <- switch(useInfo, none = NULL, pval = .pval, def = .def, 
        counts = .counts, np = paste(.def, .pval, sep = "\\\n"), 
        all = paste(.def, .pval, .counts, sep = "\\\n"))
    if (sigForAll) 
        sigNodes <- termsP.value[nodes(dag)]
    else sigNodes <- sigTerms
    if (is.null(wantedNodes)) 
        wantedNodes <- names(sigTerms)
    complete.dag <- plotFunction(dag, sigNodes = sigNodes, genNodes = names(sigTerms), 
        wantedNodes = wantedNodes, showEdges = showEdges, useFullNames = useFullNames, 
        oldSigNodes = oldSigNodes, nodeInfo = nodeInfo, sigFontSize = sigFontSize , sigBoxSize = sigBoxSize, elseFontSize = elseFontSize)
    if (swPlot && !is.null(complete.dag)) 
        plot(complete.dag)
    return(list(dag = dag, complete.dag = complete.dag))
}


environment(showSigOfNodes) <- asNamespace("topGO")
#################
## Custom plot ##
#################
GOplot_custom <- function (dag, sigNodes, dag.name = "GO terms", edgeTypes = T, 
    nodeShape.type = c("box", "circle", "ellipse", "plaintext")[3], genNodes = NULL, wantedNodes = NULL, showEdges = T, useFullNames = F, 
    oldSigNodes = NULL, nodeInfo = NULL, sigFontSize , sigBoxSize, elseFontSize) 
{
    if (!missing(sigNodes)) 
        sigNodeInd = TRUE
    else sigNodeInd = FALSE
    graphAttrs <- getDefaultAttrs(layoutType = "twopi")
    graphAttrs$cluster <- NULL
    graphAttrs$node$shape <- "ellipse"
    graphAttrs$node$fontsize <- elseFontSize
    graphAttrs$node$fixedsize <- FALSE
#    graphAttrs$node$fontcolor <- "black"
    graphAttrs$node$height <- 6
    graphAttrs$node$witdh <- 9
    nodeAttrs <- list()
    edgeAttrs <- list()
    if (is.null(nodeInfo)) {
        nodeInfo <- character(numNodes(dag))
        names(nodeInfo) <- nodes(dag)
    }
    else nodeInfo <- paste("\\\n", nodeInfo, sep = "")
    node.names <- nodes(dag)
    if (!useFullNames) 
        nodeAttrs$label <- sapply(node.names, function(x) {
            return(paste(substr(x, 4, nchar(node.names[1])), 
                nodeInfo[x], sep = ""))
        })
    else {
        nodeAttrs$label <- paste(node.names, nodeInfo, sep = "")
        names(nodeAttrs$label) <- node.names
    }
    if (!is.null(wantedNodes)) {
        diffNodes <- setdiff(wantedNodes, genNodes)
        if (length(diffNodes) > 0) {
            nodeAttrs$color[diffNodes] <- rep("lightblue", .ln <- length(diffNodes))
            nodeAttrs$shape[diffNodes] <- rep("circle", .ln)
            nodeAttrs$height[diffNodes] <- rep("0.45", .ln)
        }
    }
    if (!is.null(genNodes)) {
        nodeAttrs$color[genNodes] <- rep("lightblue", .ln <- length(genNodes))
        nodeAttrs$shape[genNodes] <- rep("box", .ln)
    }
    if (sigNodeInd) {
        if (!is.null(oldSigNodes)) {
            old.logSigNodes <- log10(sort(oldSigNodes[nodes(dag)]))
            old.range <- range(old.logSigNodes)
            logSigNodes <- log10(sort(sigNodes))
            logSigNodes[logSigNodes < old.range[1]] <- old.range[1]
            logSigNodes[logSigNodes > old.range[2]] <- old.range[2]
        }
        else old.logSigNodes <- logSigNodes <- log10(sort(sigNodes))
        sigColor <- round(logSigNodes - range(logSigNodes)[1] + 
            1)
        old.sigColor <- round(old.logSigNodes - range(old.logSigNodes)[1] + 
            1)
        mm <- max(sigColor, old.sigColor)
        sigColor <- sigColor + (mm - max(sigColor))
        colorMap <- heat.colors(mm)
        nodeAttrs$fillcolor <- unlist(lapply(sigColor, function(x) return(colorMap[x])))
    }
    if (!showEdges) 
        graphAttrs$edge$color <- "white"
    else if (edgeTypes) 
        edgeAttrs$color <- ifelse(.getEdgeWeights(dag) == 0, 
            "black", "black")
nodeAttrs$fixedsize[genNodes] <- TRUE
nodeAttrs$fontsize[genNodes] <- sigFontSize
nodeAttrs$fontcolor[genNodes] <- rep("black", .ln)
nodeAttrs$height[genNodes] <- rep(15, .ln)
nodeAttrs$width[genNodes] <- sigBoxSize


    return(agopen(graph = dag, name = dag.name, attrs = graphAttrs, 
        nodeAttrs = nodeAttrs, edgeAttrs = edgeAttrs))
}



environment(GOplot_custom) <- asNamespace("topGO")
#################




