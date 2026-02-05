function [mainTree,subTree] = extract(index,parentExpr)
%EXTRACT Extract a subtree from an encoded tree expression.
%
%   [MAINTREE,SUBTREE] = EXTRACT(INDEX,PARENTEXPR)
%
%   Input args:
%   PARENTEXPR (the parent string expression)
%   INDEX (the index in PARENTEXPR of the root node of the subtree to be
%   extracted)
%
%   Output args:
%   MAINTREE (PARENTEXPR with the removed subtree replaced by '$')
%   SUBTREE  (the extracted subtree)
%
%   Copyright (c) 2009-2015 Dominic Searson
%
%   GPTIPS 2
%
%   See also PICKNODE, GETDEPTH, GETCOMPLEXITY

endpos = numel(parentExpr);

if index < 1 || index > endpos
    error('extract: index %d out of range for expr "%s".', index, parentExpr);
end

% -------------------------------------------------------------------------
%  SAFETY: ensure INDEX is on a valid node root.
%          Valid roots are:
%            - 'x' terminals (x1,x2,...)
%            - '[' ERCs  ([...])
%            - '?' ERC token
%            - function name letter followed by '('
% -------------------------------------------------------------------------
    function tf = isValidRoot(pos)
        if pos < 1 || pos > endpos
            tf = false; return;
        end
        ch = parentExpr(pos);
        if ch == 'x' || ch == '[' || ch == '?'
            tf = true;
        elseif isletter(ch) && pos < endpos && parentExpr(pos+1) == '('
            tf = true;
        else
            tf = false;
        end
    end

if ~isValidRoot(index)
    % Snap left to the nearest valid node root (e.g. '[' for an ERC)
    found = false;
    for k = index:-1:1
        if isValidRoot(k)
            index = k;
            found = true;
            break;
        end
    end
    if ~found
        error('Malformed tree expression in EXTRACT.\nparentExpr: "%s"\nindex: %d', ...
            parentExpr, index);
    end
end

cnode  = parentExpr(index);
iplus  = index + 1;
iminus = index - 1;

if cnode == 'x'  %extracting an input terminal (x1, x2 etc.)
    
    section      = parentExpr(iplus:endpos);
    inp_comma_ind = strfind(section,',');
    inp_brack_ind = strfind(section,')');
    
    %if none found then string must consist of single input
    if isempty(inp_brack_ind) && isempty(inp_comma_ind)
        mainTree = '$';
        subTree  = parentExpr;
    else
        inp_ind   = sort([inp_brack_ind inp_comma_ind]);
        final_ind = inp_ind(1) + index;
        subTree   = parentExpr(index:final_ind-1);
        mainTree  = [parentExpr(1:iminus) '$' parentExpr(final_ind:endpos)];
    end
    
    return
    
elseif cnode == '[' %ERC
    
    cl_sbr    = strfind(parentExpr(iplus:endpos),']');
    if isempty(cl_sbr)
        error('extract: missing closing "]" for ERC in "%s".', parentExpr);
    end
    final_ind = cl_sbr(1)+index;
    subTree   = parentExpr(index:final_ind);
    mainTree  = [parentExpr(1:iminus) '$' parentExpr(final_ind+1:endpos)];
    return
    
elseif cnode=='?' %ERC token
    subTree          = cnode;
    mainTree         = parentExpr;
    mainTree(index)  = '$';
    return
    
else % otherwise extract a tree with a function node as root
    
    % subtree defined when number open brackets = number of closed brackets
    search_seg = parentExpr(index:endpos);
    
    % there must be an opening '(' after the function name
    firstOpen = strfind(search_seg,'(');
    if isempty(firstOpen)
        error('extract: malformed function node in "%s" at index %d (no opening "(").',...
              parentExpr, index);
    end
    
    firstOpen = firstOpen(1);
    balance   = 0;
    j         = 0;
    
    % walk from the first '(' onwards, find where balance returns to zero
    for t = firstOpen:numel(search_seg)
        ch = search_seg(t);
        if ch == '('
            balance = balance + 1;
        elseif ch == ')'
            balance = balance - 1;
            if balance == 0
                j = t;
                break;
            end
        end
    end
    
    if j == 0
        error('extract: mismatched parentheses in "%s" at index %d.',...
              parentExpr, index);
    end
    
    subTree  = search_seg(1:j);
    mainTree = [parentExpr(1:iminus) '$'  parentExpr(j+index:endpos)];
    return
end
end