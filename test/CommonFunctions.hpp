#ifndef COMMONFUNCTIONS_
#define COMMONFUNCTIONS_
/**
 * Functions commonly used for the simulations 
 * 
 */

void zeroFill(std::stringstream& output, unsigned value, unsigned digits)
{
    for ( unsigned d = (digits-1); d > 0; d-- )
    {
        if ( value/std::pow(10,d) < 1 )
        {
            output << "0";
        }
    }
    output << value;
}

unsigned getSeedAndAppendToFolder(std::stringstream& output, unsigned digits, unsigned multiplier = 1)
{
    unsigned seed = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-seed");
    zeroFill(output,seed,digits);
    return seed;
}

#endif
