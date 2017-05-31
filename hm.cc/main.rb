# Nazev souboru, ve kterem je chronologicky
# serazeny seznam udalosti pro simulaci.
NAZEV_SOUBORU_UDALOSTI = "udalosti.txt"

# Cena za jeden kus materialu (ocelova konstrukce luzka).
CENA_KONSTRUKCE = 10
# Cena za jeden kus materialu (matrace pro luzko).
CENA_MATRACE = 5
# Cena za jeden kus materialu (elektroinstalace pro luzko).
CENA_ELEKTROINSTALACE = 15
# Pokuta za jeden den, kdy nakladni vozidlo cekana na vylozeni.
POKUTA_PRODLENI_MATERIAL = 2
# Pokuta za jeden den, kdy nebyla objednavka kompletne vyrizena.
POKUTA_PRODLENI_OBJEDNAVKA = 5
# Prodejni cena jednoho luzka.
CENA_LUZKO = 50

class Firma
  attr_reader :rozpocet, :sklad, :den, :objednavky_vyrizene, :objednavky_k_vyrizeni

  def initialize(sklad,rozpocet,vyrobni_kapacita)
    # Kolik luzek muze firma vyrobit za jeden den.
    @vyrobni_kapacita = vyrobni_kapacita
    # Vychozi rozpocet firmy. Firma muze se muze behem simulace dostat i do zaporneho rozpoctu.
    # To simulaci nijak nevadi a neprovadejte zadne zvlastni kroky.
    @rozpocet = rozpocet
    @sklad = sklad

    # Citac dni simulace, ktery je po kazdem dni simulace zvysen o jednicku.
    @den = 0
    # Citac vyrobnich cisel luzek. Po kazdem vyrobenem luzku se zvysi o jednicku.
    @citac_vyrobnich_cisel = 1
    # Chronologicky serazeny seznam objednavek, ktere maji byt vyrizeny.
    @objednavky_k_vyrizeni = []
    # Chronologicky serazeny seznam vyrizenych objednavek.
    @objednavky_vyrizene = []
    # Uchovava informace o cekajich nakladacich s jednotlivymi druhy materialu.
    # Pod klici :konstrukce, :matrace a :elektroinstalace jsou ulozena pole
    # obsahujici pocty kusu materialu z jednotlivych nakladaku, ktere jeste nebyly
    # vylozeny, protoze nebyl dostatek mista na sklade.
    @cekajici_nakladaky = {:konstrukce => [], :matrace => [], :elektroinstalace => []}
  end

  # Nejdrive se zaplati za material, ktery byl privezen.
  # Pote se ulozi material do skladu a pripadne se zaeviduje,
  # ze nektere nakladaky se nepodarilo plne vylozit a cekaji na
  # vylozeni.

  def prijmout_material(material)
    unless material[:konstrukce].empty?
      @rozpocet -= material[:konstrukce].inject{|sum,m| sum + m } * CENA_KONSTRUKCE
    end
    unless material[:elektroinstalace].empty?
      @rozpocet -= material[:elektroinstalace].inject{|sum,m| sum + m } * CENA_ELEKTROINSTALACE
    end
    unless material[:matrace].empty?
      @rozpocet -= material[:matrace].inject{|sum,m| sum + m } * CENA_MATRACE
    end

    # IMPLEMENTUJE!
    # Prijmete do skladu jednotlive typy materialu.
    # Pokud se nepodari kompletne vylozit nakladak s materialem,
    # pak jej zaevidujte do cekajicich nakladaku dle typu materialu,
    # ktery dovezl. Pokud jiz nejake nakladaky cekaji, pak se nejdrive
    # vylozi material z nich a az po nich se pripadne vykladaji nakladaky,
    # ktere nove prijely.

    if !material[:konstrukce].empty? or !@cekajici_nakladaky[:konstrukce].empty?
      naklad =  @cekajici_nakladaky[:konstrukce] + material[:konstrukce]
      c = 0
      while c == 0 and !naklad.empty?
        pocet = naklad.shift
        odobralo = @sklad.prijmout_konstrukce(pocet)
        c = pocet - odobralo
      end
      if c != 0
        naklad.unshift(c)
      end
      @cekajici_nakladaky[:konstrukce] = naklad
    end

    if !material[:elektroinstalace].empty? or !@cekajici_nakladaky[:elektroinstalace].empty?
      naklad =  @cekajici_nakladaky[:elektroinstalace] + material[:elektroinstalace]
      c = 0
      while c == 0 and !naklad.empty?
        pocet = naklad.shift
        odobralo = @sklad.prijmout_elektroinstalace(pocet)
        c = pocet - odobralo
      end
      if c != 0
        naklad.unshift(c)
      end
      @cekajici_nakladaky[:elektroinstalace] = naklad
    end

    if !material[:matrace].empty? or !@cekajici_nakladaky[:matrace].empty?
      naklad =  @cekajici_nakladaky[:matrace] + material[:matrace]
      c = 0
      while c == 0 and !naklad.empty?
        pocet = naklad.shift
        odobralo = @sklad.prijmout_matrace(pocet)
        c = pocet - odobralo
      end
      if c != 0
        naklad.unshift(c)
      end
      @cekajici_nakladaky[:matrace] = naklad
    end

  end

  # Ulozi objednavky k vyrizeni v poradi v jakem byly
  # metode predany v parametru objednavky, ktery je typu pole.
  def prijmout_objednavky(objednavky)
    objednavky.each do |obj|
      obj[:dodana_luzka] = []
      @objednavky_k_vyrizeni << obj
    end
  end

  # IMPLEMENTUJTE!
  # Nejdrive meoda zjisti kolik luzek lze vyrobit.
  # To se urci na zaklade nasledujiciho:
  #   Nelze vyrobit vice luzek, nez se jich vejde do skladu.
  #   Nelzte vyrobit vice luzek, nez je vyrobni kapacita firmy.
  #   Pro vyrobu jednoho luzka je zapotrebi 1ks od kazdeho typu materilu.
  #     Je tedy potreba 1ks ocelove konstrukce, 1ks elektroinstalace 1ks matrace.
  # Pak probiha samotna vyroba. Tedy vytvareni instanci tridy Luzko, s tim, ze
  # nazev luzka bude LUZKO_X, kde X je hodnota z @citac_vyrobnich_cisel, ktery po
  # vyrobeni jednoho luzka zvyste o jednicku.
  # Vyrobena luzka ulozte na sklad a odeberte spotrebovany material ze skladu.
  def vyrobit_luzka
    # do promenne kolik ulozime pocet luzek k vyrobe dle uvedenych kryterii nahoru
    kolik = @sklad.volne_misto_luzka
    if @vyrobni_kapacita < kolik
      kolik = @vyrobni_kapacita
    end

    if @sklad.pocet_konstrukci < kolik
      kolik = @sklad.pocet_konstrukci
    end
    if @sklad.pocet_elektroinstalaci < kolik
      kolik = @sklad.pocet_elektroinstalaci
    end
    if @sklad.pocet_matraci < kolik
      kolik = @sklad.pocet_matraci
    end
    postilky = []
    kolik.times do
      postilky << Luzko.new(@citac_vyrobnich_cisel, @den)
      @citac_vyrobnich_cisel += 1
    end
    @sklad.prijmout_luzka(postilky)
    @sklad.odebrat_konstrukce(kolik)
    @sklad.odebrat_elektroinstalace(kolik)
    @sklad.odebrat_matrace(kolik)

  end

  # Pokud mame nejake objednavky k vyrizeni a zaroven mame na sklade pripravena luzka,
  # pak zpracovavame objednavky. Pokud je objednavka vyrizena, pak je presunuta do
  # vyrizenych objednavek.
  def zpracovat_objednavky
    while (!@objednavky_k_vyrizeni.empty? and !@sklad.luzka.empty?)
      prodano = @objednavky_k_vyrizeni[0][:pocet_ks] > @sklad.luzka.length ? @sklad.luzka.length : @objednavky_k_vyrizeni[0][:pocet_ks]
      @sklad.vyskladnit_luzka(prodano).each do |luzko|
        @objednavky_k_vyrizeni[0][:dodana_luzka] << luzko
      end
      @objednavky_k_vyrizeni[0][:pocet_ks] -= prodano
      @rozpocet += prodano * CENA_LUZKO
      if @objednavky_k_vyrizeni[0][:pocet_ks] == 0
        @objednavky_vyrizene << @objednavky_k_vyrizeni.shift
      end
    end
  end

  # IMPLEMENTUJTE
  # Za kazdou objednavku, ktera je k vyrizeni a za
  # kazdy jednotlivy cekajici nakladak odectete z firemniho rozpoctu
  # prislusnou pokutu.
  def zaplatit_pokuty
    @rozpocet -= (@cekajici_nakladaky[:konstrukce].length + @cekajici_nakladaky[:elektroinstalace].length + @cekajici_nakladaky[:matrace].length)*POKUTA_PRODLENI_MATERIAL
    @objednavky_k_vyrizeni.each{|obj| @rozpocet -= obj[:pocet_ks] * POKUTA_PRODLENI_OBJEDNAVKA}
  end

  def simuluj_den(material,objednavky)
    prijmout_material(material)
    prijmout_objednavky(objednavky)
    vyrobit_luzka
    zpracovat_objednavky
    zaplatit_pokuty
    @den += 1
  end

end

class Sklad
  attr_reader :konstrukce, :elektroinstalace, :matrace, :luzka

  def initialize(kap_konstrukce,kap_elektroinstalace,kap_matrace,kap_luzka)
    # Kapacity kolik kusu konstrukci, elektroinstalaci, matraci a luzek se vejde na sklad
    @kap_konstrukce = kap_konstrukce
    @kap_elektroinstalace = kap_elektroinstalace
    @kap_matrace = kap_matrace
    @kap_luzka = kap_luzka

    # Pocet konstrukci, elektroinstalaci a matraci ulozen na sklade
    @konstrukce = 0
    @elektroinstalace = 0
    @matrace = 0
    # Uskladnena luzka
    @luzka = []
  end

  # IMPLEMENTUJTE!
  # Odebere dany pocet kusu konstrukci ze skladu.
  # Vrati pocet odebranych kusu ze skladu.
  # Pokud bude pocet konstrukci na sklade napr. 5 a
  # parametr pocet bude 7, pak bude odebrano pouze 5
  # konstrukci.
  # Metoda vrati skutecne odebrany pocet kusu ze skladu.
  def odebrat_konstrukce(pocet)
    if pocet > @konstrukce
      n = @konstrukce
      @konstrukce = 0
      return n
    else
      @konstrukce -= pocet
      return pocet
    end
  end


  # IMPLEMENTUJTE!
  # Logika je stejna jako u odebrat_konstrukce, ale odebira
  # material typu elektroinstalace.
  def odebrat_elektroinstalace(pocet)
    if pocet > @elektroinstalace
      n = @elektroinstalace
      @elektroinstalace = 0
      return n
    else
      @elektroinstalace -= pocet
      return pocet
    end
  end


  # IMPLEMENTUJTE!
  # Logika je stejna jako u odebrat_konstrukce, ale odebira
  # material typu matrace.
  def odebrat_matrace(pocet)
    if pocet > @matrace
      n = @matrace
      @matrace = 0
      return n
    else
      @matrace -= pocet
      return pocet
    end
  end

  def pocet_konstrukci()
    return @konstrukce
  end

  def pocet_elektroinstalaci()
    return @elektroinstalace
  end

  def pocet_matraci()
    return @matrace
  end


  # IMPLEMENTUJTE!
  # Metoda prida do skladu dany pocet konstrukci,
  # ale maximalne tolik kolik se jich na sklad jeste vejde.
  # Pokud bude napr. pocet = 10 a do skladu se uz vejdou
  # pouze 3 konstrukce, pak budou prijaty pouze 3.
  # Metoda vraci pocet konstrukci, ktere se podarilo prijmout
  # do skladu.
  def prijmout_konstrukce(pocet)
    if @konstrukce+pocet > @kap_konstrukce
      k = @kap_konstrukce - @konstrukce
      @konstrukce = @kap_konstrukce
      return k
    else
      @konstrukce += pocet
      return pocet
    end
  end
  # IMPLEMENUJTE!
  # Logika je stejna jako v metode prijmout_konstrukce,
  # ale prijima elektroinstalace.
  def prijmout_elektroinstalace(pocet)
    if @elektroinstalace+pocet > @kap_elektroinstalace
      k = @kap_elektroinstalace - @elektroinstalace
      @elektroinstalace = @kap_elektroinstalace
      return k
    else
      @elektroinstalace += pocet
      return pocet
    end
  end
  # IMPLEMENUJTE!
  # Logika je stejna jako v metode prijmout_konstrukce,
  # ale prijima matrace.
  def prijmout_matrace(pocet)
    if @matrace+pocet > @kap_matrace
      k = @kap_matrace - @matrace
      @matrace = @kap_matrace
      return k
    else
      @matrace += pocet
      return pocet
    end
  end

  # Vraci pocet volnych mist pro uskladneni hotovych luzek.
  def volne_misto_luzka
    return @kap_luzka - @luzka.length
  end

  # Metoda prijme na sklad luzka, ale maximalne jich
  # prijme tolik, kolik se na sklad vejde.
  # Metoda vrati pole obsahujici NEprijata luzka.
  # Pokud jsou prijata vsechny, pak vraci nil.
  def prijmout_luzka(luzka)
    if @kap_luzka - @luzka.length >= luzka.length
      luzka.each do |l|
        @luzka << l
      end
      return nil
    else
      (@kap_luzka - @luzka.length).times do
        @luzka << luzka.shift
      end
      return luzka
    end
  end

  # IMPLEMENTUJTE!
  # Metoda vrati pole obsahujici dany pocet luzek ze skladu.
  # Luzka jsou ze skladu brana v poradi v jakem byla
  # na sklad ulozena. Tzn. prvni uskladnene luzko bude
  # vyskladneno jako prvni.
  def vyskladnit_luzka(pocet)
    pole_luzek = []
    if pocet >= @luzka.length
      pole_luzek = @luzka
      @luzka = []
      return pole_luzek
    end
    pocet.times do
      pole_luzek << @luzka.shift
    end
    return  pole_luzek
  end
end


class Luzko
  attr_reader :vyrobni_cislo,:den_vyroby

  def initialize(vyrobni_cislo,den_vyroby)
    @vyrobni_cislo = vyrobni_cislo
    @den_vyroby = den_vyroby
  end
end

sklad = Sklad.new(30,50,30,100)
firma = Firma.new(sklad,1000,30)

prijem_matrialu = {:konstrukce => [], :elektroinstalace => [], :matrace => []}
objednavky = []

File.open(NAZEV_SOUBORU_UDALOSTI,"r") do |soubor|
  # IMPLEMENTUJTE!
  # Ctete data ze vstupniho souboru a pro kazdy den simulace zavolejte na
  # firme metodu simuluj_den, ktera ocekava dva parametry: prijem_materialu a objednavky
  #
  # Prijem materilu ocekava jako asociativni pole kde pod klici :konstrukce, :elektroinstalace a :matrace
  # jsou pole obycejna obsahujici cela cisla, ktera reprezentuji jednotlive nakladaky s danym poctem kusu prislusneho
  # materialu. Pokud v dany den neprijede nakladak s danym typem materialu, pak je pod prislusnym klicem ulozeno prazdne pole.
  #
  # Objednavky ocekava jako pole obsahujici objednavky, kde jednotlive objednavky jsou reprezentovany
  # asociativnim polem, kde je pod klicem :kod ulozen kod objednavky a pod klicem :pocet_ks objednany pocet kusu luzek.
  # Pokud v dany den neni zadna objednavka, pak je toto pole prazdne.
  #
  # PRIKLAD:
  # Pro druhy den simulace (na zaklade dodanych vstupnich dat) byla metoda simuluj_den zavolana s hodnotami parametru:
  # prijem_materialu = {:konstrukce => [10,15,100], :elektroinstalace => [50], :matrace => [30]}
  # objednavky = [{:kod => "OBJEDNAVKA_3", :pocet_ks => 20},{:kod => "OBJEDNAVKA_4", :pocet_ks => 30}]

  soubor.each_line do |line|
    c = line.split(',')

    if c[0] == "DEN_SIMULACE\n"
      firma.simuluj_den(prijem_matrialu,objednavky)

      prijem_matrialu[:konstrukce] = []
      prijem_matrialu[:elektroinstalace] = []
      prijem_matrialu[:matrace] = []
      objednavky = []
    elsif c[0] == "PRIJEM_MATERIALU"
      prijem_matrialu[c[1].downcase.to_sym] << c[2].to_i
    else
      objednavky << {:kod => c[0], :pocet_ks => c[1].to_i}
    end
  end
  firma.simuluj_den(prijem_matrialu,objednavky)
end

puts "Konec simulace"
puts "Pocet dni simulace: #{firma.den}"
puts "Rozpocet: #{firma.rozpocet}"
puts "Nevyrizene objednavky:"
firma.objednavky_k_vyrizeni.each do |obj|
  puts "\t#{obj[:kod]} #{obj[:pocet_ks]}"
end
puts "Vyrizene objednavky:"
firma.objednavky_vyrizene.each do |obj|
  puts "\t#{obj[:kod]}"
  puts "\tDodana luzka:"
  obj[:dodana_luzka].each do |luzko|
    puts "\t\t#{luzko.vyrobni_cislo}, Vyrobeno: #{luzko.den_vyroby}. den"
  end
end