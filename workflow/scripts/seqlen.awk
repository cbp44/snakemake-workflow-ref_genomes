/^>/ {
    seqtotal+=seqlen
    seqlen=0
    seq+=1
    next
}
{
    seqlen += length($0)
}
END{
    print seqtotal+seqlen
}
