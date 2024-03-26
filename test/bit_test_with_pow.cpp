#include <gtest/gtest.h>

#include "../big_int/big_int.h"
#include "../big_int/big_int_ops.h"

std::pair<bint_t, bint_t> get_test_pow_values() {
    bint_t a(3ll);
    bint_t b(5ll);
    big_int::fast_pow_inplace(a, 10000);
    big_int::fast_pow_inplace(b, 5000);
    return std::make_pair(a, b);
}

TEST(BigTestWithPow, TestFirst) {
    auto [a, _] = get_test_pow_values();
    EXPECT_STREQ(a.to_string().c_str(),
        "16313501853426258743032567291811547168121324535825379939348203261918257308143190787480155630847848309673252"
        "04522323579543340558299917720385238147914536811250145319235516622439102542362884355668655965964501201417744"
        "82755299903732744254464257512355373418673876078136199372256168728620165048055931740599095204616685006631189"
        "26911571773452255850626968526251879139867085080472539640933730243410152186914328917354576854457274195562218"
        "01333774562850247067305942699911420254077317598819984248727618368529938892782529678644025299944478569418367"
        "53235217044321957858062701233883829317701989908413008615069961089447820650151634103448949458093376891568076"
        "86673462563038164792190665340124344133980763205594364754963451564072340502606377790585114123814919001637177"
        "03445738501993906023292519447111423589297856532241562834414218484289208346622787576050127600980153070303752"
        "58391578938757411924977053004696910624543699267959754563402367777343546671390726015749698343127696535571843"
        "96147587071260443947944862235744459711204473062937764153770030210332183635531818173456618022745975055313212"
        "59851442958754554729653460959719483603654687049177192762521435295750345494840363582234572877488517580950015"
        "84518373894137980953297119930921014174284067743261264500054678887365462549486586024844945359388886565427469"
        "77424368385335496083164921318601934977025095780370104307980276356857350349205866078371806065542393536101673"
        "40201798095159894698066433039150584580367424834887807101041291866733582384989962348621505030405257778984851"
        "24102638348117192369493114234118235853164050853061649366711374569853942856773247717750460509708655208935961"
        "51687017153855755197348199659070192954771308347627111052471134476325986362838585959552209645382089055182871"
        "85486674463373753321752488011840178759509406085571701014408713649553241854424148943708007471615840489591413"
        "64518020324467079610587576333456916967432938696237454108700518515906728593470612125734465720450884654606168"
        "26082579731686004585218284333452396157730036306379421822435818001505905203918209206969662326706952623512427"
        "38024046878411453510149673398340124021984004895673368930962032161379375715672756246165193339754026679596386"
        "59215909133220605726733498492533033978742423819607753371827300377836987087487817384197476988803216011863105"
        "06332869704931303076839444790968339306301273371014087248060946851793697973114432706759288546077622831002526"
        "80055484969686771028094594660366959379735464213662223119269502732122951191295294032087976312315176055595949"
        "69611631414556882788429495872883991002736918800187741475688926501861520653352191130725824176996169019955302"
        "49937735219099786758954892534365835235843156112799728164123461219817343904782402517111603206575330527850752"
        "56464299531806498590081555797994588593112435130325281125525429579708228194665879870597907749246984964418316"
        "65859508449531647268961461682978081783984704515613205261805423108407448431074693689597077268366084718170605"
        "98771730170755446473440774031371227437651048421606224757527085958515947273151027400662948161111284777828103"
        "53149948891367280078316788805117715542728510386173665806940479769590075882046523867397088266016228510759922"
        "14187436570068725378426778837088075158503976918124338805617726523648472970195080258489648338832251656689869"
        "35081274596293983121864046277268590401580209059988500511262470167150495261908136688693861324081559046336288"
        "96303709031203352240072236088249492818280907540691431995704492750442079727811783767743144697908575643299075"
        "35825881024402406110390845164010899488684333537484441046397340745191650676329414193479856244355673420728159"
        "10754484123812917487312938280670403228188813003978384081332242484646571417574404852962675165616101527367425"
        "65486950871200178839384617178045745596304576494356596488751839648129615990247199673550885429296453679677940"
        "43772309657233616251820307982977347858546060603234190916467111386784909288401074499234568347637631142260007"
        "70316931243666699425694828181155048843161380832067845480569758457751090640996007242018255400627276908188082"
        "60179552016705470132780236698974708283548110554387844688989623069609188164354747615499857401590739605947868"
        "49785741804867989184386431646185413516892583790423264876694797333847129967542517038080378286365996544477277"
        "95924596382283226723503386540591321268603222892807562509801015765174359627788357881606366119032951829868274"
        "61753994692122133028425702705865316229248268667927526676400988198559064853454493922429668979119535578320596"
        "84924226362776567353384882991042380602892093906544673162915912197128660526613470268552612893812368810630682"
        "19249064767086495184176816629077103667131505064964190910450196502178972477361881300608688593782509793781457"
        "17039689749690886189303463489571511711460151465438134713909234583347222649365693099604501635580816298496520"
        "3661519182202145414866559662218796964329217241498105206552200001"
    );
}

TEST(BigTestWithPow, TestSecond) {
    auto [_, b] = get_test_pow_values();
    EXPECT_STREQ(b.to_string().c_str(),
        "70798112610481728923856151586940575529475485103394313587298302235463672591897885235693661717820941949087897"
        "34525639653391838207739880314591126095835740247724146400971957397335974932582925939115335047391295622682040"
        "71543557621957154167022936179401953583027498888113205291131968782470891300838597837767141120435864584992782"
        "49897769718815736187243082782329229640147176436028124753382245531341376668206055393242057451777630141405462"
        "67707976418358497974514287516444661811837037264124271570952304168887253850407036602265921533670599784086381"
        "48558513504111695621640692451932474392591862086414819185411626773791894994406168428507759246036924025777055"
        "35763489639216902854903667410391791039012312202623731811269601786346292030198882360040657514640664022004942"
        "20918963695905981714074669647281359046848557531178324898469936047599147830640788861400899647898976560501012"
        "27066192624846386331984057164665045089924494201539625270697001187881367962390678914347793862826506105520522"
        "64431244602553114455145810118635386336968034465563893665924346123704717311399732341924041230923448473309935"
        "27108749236120212674650437119700471648028142556069118947380408414577374627699263517659214972183507187848497"
        "52154228451165042072134597327525335997948231177086940097213446082140059915107029942268346321886596896076937"
        "50892719901256089422722986117083506319327556202165510153623633039239256028426983308085006922536414574404177"
        "47191786881073296844297436225010128568142121796219564835210952952560722464338046736410981694424217560726644"
        "12572756068078247510287091884049267578536219177944311892558255695997942118063005962081586901411245274525575"
        "80651613972089145329366007729168154486971813631796038821675914449870059395303402413913573510312400650792013"
        "04115630512907398759840487390204883260685292053507819776468150460554434869812913994605290945663057945256291"
        "46876693677088044738908686804347196497790461267195461802633892093576728372128157045527262052597588192804596"
        "93469396531370351421002505416363522994146622818210486957034651507816836586938188210996794559074119759244226"
        "34512363637502580927278740239795409601142728413238210525169414361342477015171699998147698395327648397766260"
        "52019008163918462618920729243597086790637140819874056907380954653357507380280101201406554851578811603612887"
        "34618481989241030979080767695839374941360384352760915776597091826412000143182996828871112931417101632360487"
        "21041917157520720560769698750143353711433114361214279426411767838198064957887934495861074624324578289358208"
        "74523018002326573371792275015803915053481566654356918630304140682380843674161521948022693397290784975190390"
        "78127586705237465994641152651670971674145750993023198739621849207871768412937278551624342119496253307258773"
        "00836795578090729208678992626954221157458493787020418738362084845881383534042002828145833218292733577827993"
        "51140151928133046956866681416868222010377328757665786112914048811034184468356702964021329974343241021330927"
        "03934175161932187420211567434675744115737440152535791816634445330041968016760186511527373590612364702374070"
        "03583893248135087484990740505872980896430544530788892403011363111919570192203901996465078313817343099408643"
        "11750332699027888526776401931414364779827948392399653273629534302153643623085259212876401132505810189856324"
        "36122954142058189239059498066538905120141530496665989273342185802387118067323924396674668273908226067589036"
        "46463506527463561731104345119704946880178733563014434188375179278675060140657093731927759286933742805639023"
        "20031056203237259164457017960626894481634963085525669157505035400390625"
    );
}

TEST(BigTestWithPow, TestDivision) {
    auto [a, b] = get_test_pow_values();
    auto c = a / b;
    EXPECT_STREQ(c.to_string().c_str(),
        "23042283546710013627962331727724155339844370646143903320198607566373949023780835482979215223855990957750351"
        "66652057995501000918051976333806021936758941473695675816595831995003705157008568934976970126579157222789344"
        "71554705215314819780406582946707327996744233621218203754935783673518733882482037296914986615661160610580387"
        "86303379322545140128753724953771323169061488177583075212165885614838752217095739315907892260104901573251511"
        "98027409568618168583241887757363263339893339601157320187784602339326507646931271468576072883288529995197361"
        "35535001832489655779523201512176556234747835561084991811625346862331258860478992407341233739523413611789896"
        "37631864636672319036443272036255494675667893052281284471435406467620682916991751997823129465283886804735765"
        "62181590726352676280483931149651190747771688810899455554067410846962576053851648736183310940453715783985161"
        "33761234126433938384177044661043515545036020846161966474626959666741023657779549160971289901971186718918446"
        "61564077210783848056837726781502566813416888577041410290099728737256839027602007283317054655084551849874095"
        "00145807218366270082009517876786362737418010332149994408455931989005955593484615895032745339552923540970554"
        "1669166520976073664713516503795067110957753101859132319158077244226900725453500752494845612587246639"
    );
}